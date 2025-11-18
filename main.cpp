#include <iostream>
#include <memory>
#include <chrono>

#include <map>
#include <set>
#include <array>
#include <utility>
#include <string>
#include <stdexcept>
#include <random>
#include <algorithm>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <easy3d/renderer/drawable.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/viewer/viewer.h>
#include "easy3d/core/point_cloud.h"
#include "easy3d/renderer/drawable_points.h"
#include "easy3d/renderer/renderer.h"

#include <queue>
#include <future>
#include <ranges>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include "thread_pool.hpp"
#include "3rd_party/nlohmann/json.hpp"

namespace py = pybind11;
namespace spd = spdlog;

constexpr bool enable_debug_logging = false;

std::map<std::string, easy3d::vec3> colors = {
    {"Al", easy3d::vec3(1.0f, 0.0f, 0.0f)},
    {"N", easy3d::vec3(0.0f, 1.0f, 0.0f)},
    {"Sc", easy3d::vec3(0.0f, 0.0f, 1.0f)},
    {"O", easy3d::vec3(0.6f, 0.6f, 1.0f)},
    {"C", easy3d::vec3(1.0f, 0.6f, 0.3f)},
    {"I", easy3d::vec3(0.8f, 0.3f, 1.0f)},
    {"Mo", easy3d::vec3(1.0f, 0.7f, 0.2f)},
};


class CellViewer {
    std::vector<easy3d::vec3> m_atomic_positions;
    std::map<int, std::pair<std::string, std::array<float, 3>>> *m_ptr_structure_dict;
    std::map<std::string, easy3d::vec3> m_unique_element_colors;
    easy3d::PointCloud m_atoms;

public:
    easy3d::Viewer viewer;

    CellViewer(std::map<int, std::pair<std::string, std::array<float, 3>>> *structure_dict_ptr) {
        this ->m_ptr_structure_dict = structure_dict_ptr;
    }

    void pick_colors() {
        //firstly access all element values
        std::vector<std::string> elements;
        for (auto const& [key, val] : *m_ptr_structure_dict) {
            elements.push_back(val.first);
        }
        std::set<std::string> element_set;
        for (auto & element : elements) {
            element_set.insert(element);
        }
        elements.assign( element_set.begin(), element_set.end());

        for (auto const& element : elements) {
            m_unique_element_colors[element] = colors[element];
        }
    }

    void initialize_viewer() {
        spd::default_logger()->log(spd::level::info, "Initializing viewer");
        viewer.resize(800, 800);
        auto colors = m_atoms.add_vertex_property<easy3d::vec3>("v:color");  // create vertex color property

        pick_colors();

        for (int i=0; i < m_ptr_structure_dict->size(); i++) {
            auto &entry = m_ptr_structure_dict->at(i);
            auto vertex = m_atoms.add_vertex(easy3d::vec3({entry.second[0], entry.second[1], entry.second[2]}));// z = 0: all points are on XY plane
            colors[vertex] = m_unique_element_colors[entry.first];
            m_atomic_positions.push_back({entry.second[0], entry.second[1], entry.second[2]});
        }
        viewer.add_model(&m_atoms);
        auto drawable = m_atoms.renderer()->get_points_drawable("vertices");
        drawable->set_impostor_type(easy3d::PointsDrawable::SPHERE);
        drawable->set_point_size(50);

        // Compute the bounding box.
        auto bbox_drawable = new easy3d::LinesDrawable("bbox");
        const easy3d::Box3 &box = easy3d::geom::bounding_box<easy3d::Box3, std::vector<easy3d::vec3>>(m_atomic_positions);
        float xmin = box.min_coord(0);
        float xmax = box.max_coord(1);
        float ymin = box.min_coord(0);
        float ymax = box.max_coord(1);
        float zmin = box.min_coord(0);
        float zmax = box.max_coord(1);
        // The eight vertices of the bounding box.
        const std::vector<easy3d::vec3> bbox_points = {
            easy3d::vec3(xmin, ymin, zmax), easy3d::vec3(xmax, ymin, zmax),
            easy3d::vec3(xmin, ymax, zmax), easy3d::vec3(xmax, ymax, zmax),
            easy3d::vec3(xmin, ymin, zmin), easy3d::vec3(xmax, ymin, zmin),
            easy3d::vec3(xmin, ymax, zmin), easy3d::vec3(xmax, ymax, zmin)
    };
        // The vertex indices of the twelve line segments of the bounding box.
        const std::vector<unsigned int> bbox_indices = {
            0, 1, 2, 3, 4, 5, 6, 7,
            0, 2, 4, 6, 1, 3, 5, 7,
            0, 4, 2, 6, 1, 5, 3, 7
    };
        // Upload the vertex positions of the bounding box to the GPU.
        bbox_drawable->update_vertex_buffer(bbox_points);
        // Upload the vertex indices of the bounding box to the GPU.
        bbox_drawable->update_element_buffer(bbox_indices);
        // Set a color for the edges of the bounding box (here we want a blue color).
        bbox_drawable->set_uniform_coloring(easy3d::vec4(0.0f, 0.0f, 1.0f, 1.0f));    // r, g, b, a
        // Set the width of the edges (here 5 pixels).
        bbox_drawable->set_line_width(5.0f);
        // Add the drawable to the viewer
        viewer.add_drawable(std::shared_ptr<easy3d::LinesDrawable>(bbox_drawable));

        viewer.fit_screen();
    }
    /*
    void update_colors() {
        auto vertices = m_atoms.get_vertex_property<std::vector<easy3d::vec3>>("v:color");
        for (int i = 0; i < m_ptr_structure_dict->size(); i++) {
            auto &entry = m_ptr_structure_dict->at(i);
            colors[] = m_unique_element_colors[entry.first];
    }
    */

    void open_viewer() {
        viewer.run();
        viewer.exit();
    }

    void update_color_buffer() {
        size_t index = 0;
        auto drawable = m_atoms.renderer()->get_points_drawable("vertices");
        auto colors = m_atoms.vertex_property<easy3d::vec3>("v:color");
        for (auto vertex : m_atoms.vertices()) {
            auto &entry = m_ptr_structure_dict->at(index);
            colors[vertex] = m_unique_element_colors[entry.first];
            ++index;
        }
        if constexpr (enable_debug_logging) {
            spd::debug("Updating color buffer");
        }
        drawable->set_coloring(easy3d::State::COLOR_PROPERTY, easy3d::State::VERTEX, "v:color");
        drawable->update();
    }
};

class Agent {
    int m_state = 0;

    public:
    std::vector<int> *neighbors_ptr;

    explicit Agent(std::vector<int> *neighbors_ptr) {
        this -> neighbors_ptr = neighbors_ptr;
    }

    void reset() {
        m_state = 0;//so 0 maps to init element before doping
    }
    void infect() {
        m_state = 1;
    }

    int get_status() const {
        if (m_state == 1) {
            //reset();
        }
        return m_state;
    }
};

class Grid {
    inline static thread_local std::default_random_engine m_generator;
    inline static std::uniform_int_distribution<int> m_uniform_int_distribution;
    std::map<int, std::vector<int>> *neighbor_dict_ptr;
    std::map<int, std::vector<int>> neighbor_dict_init;
    std::vector<int> *m_to_infect_ptr {};

public:
    std::map<int, std::unique_ptr<Agent>> agent_unique_ptrs;
    explicit Grid(std::map<int, std::vector<int>> *neighbor_dict_ptr) {
        // to infect has size 0
        this -> neighbor_dict_ptr = neighbor_dict_ptr;
    }

    const std::map<int, std::unique_ptr<Agent>> &get_agent_unique_ptrs() const {
        //return a const ref to the vector of unique ptrs
        return agent_unique_ptrs;
    }

    void init_simulation (std::vector<int> *to_infect_ptr) {
        neighbor_dict_init = *neighbor_dict_ptr;

        m_to_infect_ptr = to_infect_ptr;
        spd::info("Initializing simulation");
        for (int i = 0; i < neighbor_dict_ptr->size(); i++) {
            agent_unique_ptrs[i] = std::make_unique<Agent>(&neighbor_dict_ptr->at(i));
            // if agent is in to_infect vector
            if (std::find(to_infect_ptr->begin(), to_infect_ptr->end(), i) != to_infect_ptr->end()) {
                agent_unique_ptrs[i]->infect();
                //spd::default_logger()->log(spd::level::debug, std::format("Infecting agent {}", i));
            }
            else {
                agent_unique_ptrs[i]->reset();
            }
        }
    }

    void draw_neighbors_to_infect_next(std::map<int, int> *next_infected_ptr) {
        for (int i = 0; i < agent_unique_ptrs.size(); i++) {
            // check if the current agent i is susceptible to the virus
            //std::cout << agent_infected << std::endl;
            bool const agent_infected = (agent_unique_ptrs[i]->get_status() == 1);

            // reference to neighbors_ptr
            auto &neighbors_ptr = agent_unique_ptrs[i]->neighbors_ptr;
            if constexpr (enable_debug_logging) {
                spd::debug("Agent {} has {} neighbors", i, neighbors_ptr->size());
            }
            std::uniform_int_distribution<int> m_uniform_int_distribution(0, neighbors_ptr->size()-1);
            // pick a random neighbor from neighbors_ptr's vector that is not in temp
            int random_neighbor_index = 2;
            bool neighbor_not_susceptible = true;

            if (neighbors_ptr->size() != 0) {
                for (int j = 0; neighbor_not_susceptible and j<neighbor_dict_init.at(i).size(); j++) {
                    if constexpr (enable_debug_logging) {
                        spd::debug("loop");
                    }
                    random_neighbor_index = m_uniform_int_distribution(m_generator);
                    if (agent_unique_ptrs.at(neighbors_ptr->at(random_neighbor_index))->get_status() == 0) {
                        neighbor_not_susceptible = false;
                    }
                }

                if (neighbor_not_susceptible) {
                    random_neighbor_index = -1;
                }
            }
            else {
                random_neighbor_index = -1;
            }


            if constexpr (enable_debug_logging) {
                spd::debug("Picked {} as neighbor", random_neighbor_index);
            }


            //std::cout << (m_picked_swaps.end() != std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) << std::endl;

            //chose random neighbor of agent i and check if susceptible

            int random_neighbor;
            if (random_neighbor_index == -1) {
                if constexpr (enable_debug_logging) {
                    spd::debug("Received -1");
                }
                random_neighbor = i;
                if (agent_infected) {
                    (*next_infected_ptr)[i] = random_neighbor;
                }

            }
            else {
                if constexpr (enable_debug_logging) {
                    spd::debug("Received ! -1");
                }
                int random_neighbor = neighbors_ptr->at(random_neighbor_index);
                if (agent_infected) {
                    if constexpr (enable_debug_logging) {
                        spd::debug("Infected agent {} has susceptible neighbor {}", i, random_neighbor);
                    }
                    (*next_infected_ptr)[i] = random_neighbor;
                }
            }
        }
    }

    void update_agent_states(std::vector<int> *next_infected_ptr) {
        for (int i = 0; i < agent_unique_ptrs.size(); i++) {
            bool const agent_in_agents_to_infect = next_infected_ptr->end() != std::find(next_infected_ptr->begin(), next_infected_ptr->end(), i);
            bool const agent_infected = (agent_unique_ptrs[i]->get_status() == 1);
            if (agent_in_agents_to_infect) {
                //std::cout << "Agent " << i << " is infected right now." << std::endl;
                agent_unique_ptrs[i]->infect();
            }
            else if (agent_infected) {
                agent_unique_ptrs[i]->reset();
                //std::cout << "Agent " << i << " is reset right now." << std::endl;
            }
        }
    }

    void map_next_infected_to_causes(std::map<int, int> *next_infected_ptr, std::map<int, std::vector<int>> *causes_of_infection_ptr) {
        for (auto &[agent_index, neighbor_index] : *next_infected_ptr) {
            auto &vector_ref = (*causes_of_infection_ptr)[neighbor_index];
            if (vector_ref.end() == std::find(vector_ref.begin(), vector_ref.end(), agent_index)) {
                vector_ref.push_back(agent_index);
                if constexpr (enable_debug_logging) {
                    spd::debug("Added {} to have caused infection of neighbor {}, Size: {}", agent_index, neighbor_index, vector_ref.size());
                }
            }
        }
    }

    void map_next_infected_to_hits(std::map<int, int> *hits_ptr, std::map<int, std::vector<int>> *causes_of_infection_ptr) {
        for (auto &[key, val] : *causes_of_infection_ptr) {
            (*hits_ptr)[key] = val.size();
            if constexpr (enable_debug_logging) {
                spd::debug("Number of hits of agent {}: {}", key, hits_ptr->at(key));
            }
        }
    }

    void update_simulation() {
        std::vector<int> next_infected_vector = {};
        // cause_of_infection maps the next infected agent to the pointers of agent pointers that caused the next infection (one to many)
        // next_infected maps the infected agent to the next infected agent (one to one)
        std::map<int, std::vector<int>> causes_of_infection = {};
        std::map<int, int> next_infected = {};
        std::map<int, int> hits = {};

        draw_neighbors_to_infect_next(&next_infected);
        map_next_infected_to_causes(&next_infected, &causes_of_infection);
        map_next_infected_to_hits(&hits, &causes_of_infection);

        // iterate through hits and if val >1 use key to access val of causes_of_infection
        // val = vector of agent indices
        // use these to remove the val of causes_of_infection from almost all neighbor dicts

        // get boolean true if any entry in collisions is set to true
        auto hits_iterator = hits.begin();
        int it = 0;
        while (std::any_of(hits.begin(), hits.end(), [](auto &key_value_pair) {return key_value_pair.second > 1;})) {
            if constexpr (enable_debug_logging) {
                spd::debug("Iteration {}", it);
            }
            for (;hits_iterator != hits.end(); ++hits_iterator) {
                if (hits_iterator->second > 1) {
                    if constexpr (enable_debug_logging) {
                        spd::debug("Collision event at lattice site {}", hits_iterator->first);
                    }
                    for (auto i : causes_of_infection.at(hits_iterator->first)) {
                        // now remove key in agent i's neighbors vector
                        if constexpr (enable_debug_logging) {
                            spd::debug("Removing {} as agent {}'s neighbor", hits_iterator->first, i);
                        }
                        auto &ptr_to_neighbors = agent_unique_ptrs.at(i)->neighbors_ptr;
                        ptr_to_neighbors->erase(std::remove(ptr_to_neighbors->begin(), ptr_to_neighbors->end(), hits_iterator->first), ptr_to_neighbors->end());
                    }
                }
                else {}
            }

            causes_of_infection.clear();
            next_infected.clear();
            hits.clear();
            draw_neighbors_to_infect_next(&next_infected);
            map_next_infected_to_causes(&next_infected, &causes_of_infection);
            map_next_infected_to_hits(&hits, &causes_of_infection);

            hits_iterator = hits.begin();
            ++it;
        }

        // fill next_infected_vector

        for (auto &[key, val] : next_infected) {
            next_infected_vector.push_back(val);
        }
        if constexpr (enable_debug_logging) {
            spd::debug("Number of next infections: {}", next_infected_vector.size());
        }

        //agents_to_infect will be known after the collision detection

        //if (agent_infected) {
        //    agents_to_infect_next_ptr->push_back(random_neighbor);
        //    spd::info(agents_to_infect_next_ptr->size());
        //}
        update_agent_states(&next_infected_vector);
        //std::this_thread::sleep_for(std::chrono::milliseconds(1));

        // we have to reset the agent's neighbor vectors individually based on the init map
        // Reset each agent's neighbors to the initial neighbors
        for (auto &[index, agent_ptr] : agent_unique_ptrs) {
            *(agent_ptr->neighbors_ptr) = neighbor_dict_init[index];
        }
    }
};


class CellOperator {
    inline static thread_local std::default_random_engine m_generator;
    std::vector<int> m_picked_swaps;
    std::map<int, std::pair<std::string, std::array<float, 3>>> *m_ptr_structure_dict;
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict_init;
    std::map<int, std::vector<int>> *m_ptr_neighbor_dict;
    Grid m_grid;
    std::vector<std::unique_ptr<Agent>> m_agent_unique_ptrs;
    std::string m_element;
    std::vector<std::string> m_config;
    std::vector<std::vector<std::string>> m_configurations = {};
    int m_n_configurations;

    public:
    CellOperator(std::map<int, std::pair<std::string, std::array<float, 3>>> *structure_dict_ptr, std::map<int, std::vector<int>> *neighbor_dict_ptr)
        :  m_ptr_structure_dict(structure_dict_ptr), m_ptr_neighbor_dict(neighbor_dict_ptr), m_grid(m_ptr_neighbor_dict)
    {
    }
    //highly problematic, copy the dict and dope the copy
    //we need the undoped dict for the virus simulation

    void backup_dict() {
        m_structure_dict_init = *m_ptr_structure_dict;
        spd::info("Copying structure dict");
    }

    void dope_cell(int n_swaps, std::string element) {
        this->m_element = element;

        std::uniform_int_distribution<int> m_uniform_int_distribution(0, m_ptr_structure_dict->size()-1);
        if (n_swaps >= m_ptr_structure_dict->size()) {
            throw std::invalid_argument("n_swaps must not be >= than the number of atoms! Don't be lazy and build the pure cell yourself...");
            }
        spd::info("Swapping atoms");
        for (int i = 0; m_picked_swaps.size() < n_swaps; i++) {
            int temp = m_uniform_int_distribution(m_generator);
            //std::cout << (m_picked_swaps.end() != std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) << std::endl;
            if (m_picked_swaps.end() == std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) {
                m_picked_swaps.push_back(temp);
                }
            }
        for (int index : m_picked_swaps) {
            m_ptr_structure_dict->at(index).first = element;
        }
        m_grid.init_simulation(&m_picked_swaps);
    }

    void update_m_structure_dict() {
        auto &unique_ptrs = m_grid.get_agent_unique_ptrs();
        for (const auto& [key, val] : *m_ptr_structure_dict) {
            int status = unique_ptrs.at(key)->get_status();
            //std::cout << "status: " << status << std::endl;

            if (status == 0) {
                //here is a mistake!!!!
                m_ptr_structure_dict->at(key).first = m_structure_dict_init.at(key).first;
                if constexpr (enable_debug_logging) {
                    spd::debug("Agent {} -> {}", key, m_ptr_structure_dict->at(key).first);
                }
            }
            else if (status == 1) {
                m_ptr_structure_dict->at(key).first = m_element;
                if constexpr (enable_debug_logging) {
                    spd::debug("Agent {} -> {}", key, m_element);
                }
            }
            else {
                if constexpr (enable_debug_logging) {
                    spd::warn("Status for agent {} is unknown.", key);
                }
            }

        }
    }

    void get_n_configurations() {
        int n_atoms = m_ptr_structure_dict->size();
        int n_swaps = m_picked_swaps.size();
        int prod_1 = 1, prod_2 = 1;
        for (int i = 0; i < n_swaps; i++) {
            prod_1 = prod_1*(n_atoms - i);
            prod_2 = prod_2*(n_swaps - i);
        }
        m_n_configurations =  prod_1 / prod_2;
        //int const n_conf
        //spd::info("Configuration space consists of {} possible configurations.", n_conf);
    }

    void run_simulation(const std::unique_ptr<CellViewer> &ptr_cell_viewer) {
        int n_found_configs = 0;
        get_n_configurations();
        // to_indext = picked_swaps
        spd::info("Truncating after {} configurations have been found", m_n_configurations);

        int i = 0;
        while (n_found_configs < m_n_configurations) {
            if (i%1000 == 0) {
                spd::info("Simulation step {}", i);
                spd::info("Number of found configurations {}", n_found_configs);
            }
            m_grid.update_simulation();
            ptr_cell_viewer->update_color_buffer();
            ptr_cell_viewer->viewer.update();
            update_m_structure_dict();

            //here write every m_structure_dict to file or append elements to vector
            // check vector of elements against vector of vector of elements

            m_config.clear();
            for (auto &[key, val] : *m_ptr_structure_dict) {
                auto &element = val.first;
                m_config.emplace_back(element);
            }

            bool exists = std::any_of(
    m_configurations.begin(), m_configurations.end(),
    [&](const auto& config) {
        return std::equal(config.begin(), config.end(), m_config.begin());
    }
);
            if (!exists) {
                m_configurations.emplace_back(m_config);
                ++n_found_configs;
            }



            // update m_structure_dict using its ptr
            // iterate over agent states and if state=0 take element from m_structure_dict_init
            // if state=1 take dopant element
            ++i;
        }
        spd::info("Simulation complete");
        spd::info("Found {} configurations after {} simulation steps", n_found_configs, i);
        // initialize the simulation grid with to_infect and pointers to index of structure_dict or neighbor_dict
        //
    }
};

class Framework {
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict;
    std::map<int, std::vector<int>> m_neighbor_dict;

    CellOperator m_cell_operator;
    CellViewer m_cell_viewer;

    public:
    ThreadPool thread_pool{1};
    std::queue<std::future<void>> results;

    py::dict python_structure_dict;
    py::dict python_neighbor_dict;
    // pass the address of the map to the individual objects
    // the map is still owned by Framework
    Framework(const py::dict &structure_dict, const py::dict &neighbor_dict, const py::kwargs& kwargs)
        : m_cell_operator(&m_structure_dict, &m_neighbor_dict), m_cell_viewer(&m_structure_dict)
    {
        if (not kwargs.contains("crystal_system")) {
            throw std::invalid_argument("Specify keyword argument 'crystal_system'!");
        }
        if (kwargs.contains("crystal_system") and kwargs["crystal_system"].cast<std::string>() != "cubic") {
            throw std::invalid_argument("Crystal system must be a cubic!");
        }
        python_structure_dict = structure_dict;
        python_neighbor_dict = neighbor_dict;
        parse_dicts();
        m_cell_operator.backup_dict();
    }

    void parse_dicts() {
        // look for general, crystal structure = cubic, if not: break
        // parse atoms into map key: (atom_index, element), value: array[3] of fractional coords
        //map will be used later

        for (auto const& [key, val] : python_structure_dict) {
            m_structure_dict[key.cast<int>()] = {val["element"].cast<std::string>(), val["frac_coord"].cast<std::array<float, 3>>()};
        }
        for (auto const& [key, val] : python_neighbor_dict) {
            m_neighbor_dict[key.cast<int>()] = {val.cast<std::vector<int>>()};
        }
    }

    void dope_cell(int n_swaps, const std::string &element) {
        m_cell_operator.dope_cell(n_swaps, element);
    }

    void display_cell() {
        m_cell_viewer.open_viewer();
    }

    void run_simulation() {
        m_cell_viewer.initialize_viewer();
        auto m_cell_viewer_unique_ptr = std::unique_ptr<CellViewer>(&m_cell_viewer);
        results.emplace(thread_pool.AddTask([this, &m_cell_viewer_unique_ptr]() { this->m_cell_operator.run_simulation(m_cell_viewer_unique_ptr); }));
        m_cell_viewer.open_viewer();

        //m_cell_operator.run_simulation();
    }
};

PYBIND11_MODULE(finding_generators, m) {
    m.doc() = "Python bindings for functions in Python";
    py::class_<Framework>(m, "Framework")
    .def(py::init<py::dict, py::dict, py::kwargs>())
    .def("dope_cell", &Framework::dope_cell)
    .def("display_cell", &Framework::display_cell)
    .def("run_simulation", &Framework::run_simulation);
}
