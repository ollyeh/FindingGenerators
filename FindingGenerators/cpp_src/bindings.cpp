#include <iostream>
#include <memory>
#include <chrono>
#include <iomanip>

#include <map>
#include <set>
#include <array>
#include <utility>
#include <string>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <fstream>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#if GRAPHICS
#include <easy3d/renderer/drawable.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/viewer/viewer.h>
#include "easy3d/core/point_cloud.h"
#include "easy3d/renderer/drawable_points.h"
#include "easy3d/renderer/renderer.h"
#endif

#include <queue>
#include <future>
#include <ranges>

#include "spdlog/spdlog.h"

#include "thread_pool.hpp"
#include "nlohmann/json.hpp"
#include "spdlog/fmt/bundled/chrono.h"

namespace py = pybind11;
namespace spd = spdlog;

constexpr bool enable_debug_logging = false;

unsigned long long binomial(unsigned n, unsigned k) {
    if (k > n) return 0;
    if (k > n - k) k = n - k;

    __uint128_t result = 1;
    for (unsigned i = 1; i <= k; ++i) {
        result = result * (n - (k - i));
        result /= i;
    }
    return (unsigned long long) result;
}

#if GRAPHICS
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
    std::shared_ptr<easy3d::PointCloud> m_atoms_shared_ptr;
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
            m_atomic_positions.emplace_back(entry.second[0], entry.second[1], entry.second[2]);
        }
        m_atoms_shared_ptr = std::make_shared<easy3d::PointCloud>(m_atoms);
        viewer.add_model(m_atoms_shared_ptr);
        auto drawable = viewer.current_model()->renderer()->get_points_drawable("vertices");
        //auto drawable = m_atoms.renderer()->get_points_drawable("vertices");
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
    }

    void close_viewer() {
        viewer.exit();
    }

    void show_property_stats() {
        int index = 0;
        //viewer.current_model()->property_stats(std::cout);
        auto positions = m_atoms_shared_ptr->get_vertex_property<easy3d::vec3>("v:point");
        for (auto vertex : m_atoms.vertices()) {
            auto &entry = m_ptr_structure_dict->at(index);
            positions[vertex] = easy3d::vec3(entry.second[0], entry.second[1], entry.second[2]);
            ++index;
        }
        viewer.current_model()->renderer()->update();
    }

    void update_color_buffer() {
        int index = 0;
        //auto colors = m_atoms.vertex_property<easy3d::vec3>("v:color");
        auto colors = m_atoms_shared_ptr->get_vertex_property<easy3d::vec3>("v:color");
        for (auto vertex : m_atoms.vertices()) {
            auto &entry = m_ptr_structure_dict->at(index);
            colors[vertex] = m_unique_element_colors[entry.first];
            ++index;
        }
        viewer.current_model()->renderer()->update();
        if constexpr (enable_debug_logging) {
            spd::debug("Updating color buffer");
        }
    }
};
#else
class CellViewer;
#endif

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
            int random_neighbor_index;
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
        int it = 0;
        while (std::any_of(hits.begin(), hits.end(), [](auto &key_value_pair) {return key_value_pair.second > 1;})) {
            if constexpr (enable_debug_logging) {
                spd::debug("Iteration {}", it);
            }
            for (auto &[key, val] : hits) {
                if (val > 1) {
                    if constexpr (enable_debug_logging) {
                        spd::debug("Collision event at lattice site {}", key);
                    }
                    for (auto i : causes_of_infection.at(key)) {
                        // now remove key in agent i's neighbors vector
                        if constexpr (enable_debug_logging) {
                            spd::debug("Removing {} as agent {}'s neighbor", key, i);
                        }
                        auto &ptr_to_neighbors = agent_unique_ptrs.at(i)->neighbors_ptr;
                        ptr_to_neighbors->erase(std::remove(ptr_to_neighbors->begin(), ptr_to_neighbors->end(), key), ptr_to_neighbors->end());
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

class AtomPermutator{
    inline static thread_local std::default_random_engine m_generator;
    std::vector<int> m_picked_swaps;
    std::map<int, std::pair<std::string, std::array<float, 3>>> *m_ptr_structure_dict;
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict_init;
    std::string m_element;
    std::vector<int> m_config;
    int m_n_configurations{};

    public:
    std::vector<std::map<int, std::pair<std::string, std::array<float, 3>>>> configurations = {};
    nlohmann::json root;

    explicit AtomPermutator(std::map<int, std::pair<std::string, std::array<float, 3>>> *structure_dict_ptr)
        :  m_ptr_structure_dict(structure_dict_ptr)
    {
    }
    //highly problematic, copy the dict and dope the copy
    //we need the undoped dict for the virus simulation

    void backup_dict() {
        m_structure_dict_init = *m_ptr_structure_dict;
        spd::info("Copying structure dict");
    }

    void dope_cell(int n_swaps, const std::string& element) {
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
    }

    void update_m_structure_dict() {
        int i = 0;
        for (const auto& [key, val] : *m_ptr_structure_dict) {
            //std::cout << "status: " << status << std::endl;

            if (m_config[i] == 0) {
                m_ptr_structure_dict->at(key).first = m_structure_dict_init.at(key).first;
                if constexpr (enable_debug_logging) {
                    spd::debug("Agent {} -> {}", key, m_ptr_structure_dict->at(key).first);
                }
            }
            else if (m_config[i] == 1) {
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
            ++i;
        }
    }

    void run_permutation(CellViewer *ptr_cell_viewer, const std::string *path) {
        m_config.clear();

        auto &cell_viewer = *ptr_cell_viewer;

        m_n_configurations = binomial(m_ptr_structure_dict->size(), m_picked_swaps.size());
        spd::info("Considering all {} permutations", m_n_configurations);

        for (auto &[key, val] : *m_ptr_structure_dict) {
            auto &element = val.first;
            if (element == m_structure_dict_init[key].first) {
                m_config.emplace_back(0);
            }
            else {
                m_config.emplace_back(1);
            }
        }

        std::sort(m_config.begin(), m_config.end());

        int configuration_index = 0;
        do {
            //std::this_thread::sleep_for(std::chrono::milliseconds(400));

             if (configuration_index%1000 == 0) {
                spd::info("Permutation step {} / {} %", configuration_index, static_cast<float>(configuration_index)/static_cast<float>(m_n_configurations)*100);
            }

            update_m_structure_dict();
            nlohmann::json configuration;
            int atom_index = 0;
            for (auto &[key, val] : *m_ptr_structure_dict) {
                configuration[std::to_string(atom_index)] = {{"element", val.first}, {"frac_coord", val.second}};
                ++atom_index;
            }
            root[std::to_string(configuration_index)] = configuration;
            #if GRAPHICS
            cell_viewer.update_color_buffer();
            cell_viewer.viewer.update();
            #endif
            ++configuration_index;
        } while (std::next_permutation(m_config.begin(), m_config.end()));
        spd::info("Done finding {} permutations", configuration_index);
        spd::info("Dumping configurations to json");
        std::ofstream(*path) << root.dump(1);
    }
};

class AtomPermutatorFramework {
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict;
    std::map<int, std::vector<int>> m_neighbor_dict;

    AtomPermutator m_atom_permutator;
    #if GRAPHICS
    CellViewer m_cell_viewer;
    #endif

    public:
    ThreadPool thread_pool{1};
    std::queue<std::future<void>> results;

    py::dict python_structure_dict;
    py::dict python_neighbor_dict;
    py::dict python_configuration_dict = {};
    // pass the address of the map to the individual objects
    // the map is still owned by Framework
    AtomPermutatorFramework(const py::dict &structure_dict)
        : m_atom_permutator(&m_structure_dict)
    #if GRAPHICS
    , m_cell_viewer(&m_structure_dict)
    #endif
    {
        python_structure_dict = structure_dict;
        parse_dicts();
        m_atom_permutator.backup_dict();
    }

    void parse_dicts() {
        // look for general, crystal structure = cubic, if not: break
        // parse atoms into map key: (atom_index, element), value: array[3] of fractional coords
        //map will be used later

        for (auto const& [key, val] : python_structure_dict) {
            m_structure_dict[key.cast<int>()] = {val["element"].cast<std::string>(), val["frac_coord"].cast<std::array<float, 3>>()};
        }
    }

    void dope_cell(int n_swaps, const std::string &element) {
        m_atom_permutator.dope_cell(n_swaps, element);
    }

    #if GRAPHICS
    void run_permutation(const std::string &path) {
        m_cell_viewer.initialize_viewer();
        auto m_cell_viewer_ptr = &m_cell_viewer;
        results.emplace(thread_pool.AddTask([this, m_cell_viewer_ptr, &path] { this->m_atom_permutator.run_permutation(m_cell_viewer_ptr, &path); }));
        m_cell_viewer.open_viewer();
    }
    #else
    void run_permutation(const std::string &path) {
        m_atom_permutator.run_permutation(nullptr, &path);
    }
    #endif

};

class VirusSimulator {
    inline static thread_local std::default_random_engine m_generator;
    std::vector<int> m_picked_swaps;
    std::map<int, std::pair<std::string, std::array<float, 3>>> *m_ptr_structure_dict;
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict_init;
    std::map<int, std::vector<int>> *m_ptr_neighbor_dict;
    Grid m_grid;
    std::vector<std::unique_ptr<Agent>> m_agent_unique_ptrs;
    std::string m_element;
    std::vector<int> m_config;
    std::string m_config_str = "";
    std::vector<std::string> m_configurations = {};
    int m_n_configurations;

    public:
    VirusSimulator(std::map<int, std::pair<std::string, std::array<float, 3>>> *structure_dict_ptr, std::map<int, std::vector<int>> *neighbor_dict_ptr)
        :  m_ptr_structure_dict(structure_dict_ptr), m_ptr_neighbor_dict(neighbor_dict_ptr), m_grid(m_ptr_neighbor_dict)
    {
    }
    //highly problematic, copy the dict and dope the copy
    //we need the undoped dict for the virus simulation

    void backup_dict() {
        m_structure_dict_init = *m_ptr_structure_dict;
        spd::info("Copying structure dict");
    }

    void dope_cell(int n_swaps, const std::string& element) {
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

    void run_simulation(CellViewer *ptr_cell_viewer) {
        auto &cell_viewer = *ptr_cell_viewer;
        int n_found_configs = 0;
        m_n_configurations = binomial(m_ptr_structure_dict->size(), m_picked_swaps.size());
        // to_indext = picked_swaps
        spd::info("Truncating after {} configurations have been found", m_n_configurations);

        int i = 0;
        while (n_found_configs < m_n_configurations) {
            m_config_str.clear();

            if (i%1000 == 0) {
                spd::info("Simulation step {}", i);
                spd::info("Number of found configurations {}", n_found_configs);
            }
            m_grid.update_simulation();
#if GRAPHICS
            cell_viewer.update_color_buffer();
            cell_viewer.viewer.update();
#endif
            update_m_structure_dict();

            //here write every m_structure_dict to file or append elements to vector
            // check vector of elements against vector of vector of elements

            m_config.clear();

            for (auto &[key, val] : *m_ptr_structure_dict) {
                auto &element = val.first;
                m_config_str = m_config_str + element;
            }

            bool exists = std::any_of(m_configurations.begin(), m_configurations.end(), [&](const auto& config) {return std::equal(config.begin(), config.end(), m_config_str.begin());});
            if (!exists) {
                m_configurations.emplace_back(m_config_str);
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

class VirusSimulatorFramework {
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict;
    std::map<int, std::vector<int>> m_neighbor_dict;

    VirusSimulator m_virus_simulator;
#if GRAPHICS
    CellViewer m_cell_viewer;
#endif

    public:
    ThreadPool thread_pool{1};
    std::queue<std::future<void>> results;

    py::dict python_structure_dict;
    py::dict python_neighbor_dict;
    py::dict python_configuration_dict = {};
    // pass the address of the map to the individual objects
    // the map is still owned by Framework
    VirusSimulatorFramework(const py::dict &structure_dict, const py::dict &neighbor_dict)
        : m_virus_simulator(&m_structure_dict, &m_neighbor_dict)
#if GRAPHICS
    , m_cell_viewer(&m_structure_dict)
#endif
    {
        python_structure_dict = structure_dict;
        python_neighbor_dict = neighbor_dict;
        parse_dicts();
        m_virus_simulator.backup_dict();
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
        m_virus_simulator.dope_cell(n_swaps, element);
    }
#if GRAPHICS
    void run_simulation() {
        m_cell_viewer.initialize_viewer();
        auto m_cell_viewer_ptr = &m_cell_viewer;
        results.emplace(thread_pool.AddTask([this, m_cell_viewer_ptr]() { this->m_virus_simulator.run_simulation(m_cell_viewer_ptr); }));
        m_cell_viewer.open_viewer();
        //return python_configuration_dict;

        //m_cell_operator.run_simulation();
    }
#else
    void run_simulation() {
        m_virus_simulator.run_simulation(nullptr);
    }
    #endif
};

class GeneratorFinder {
    // this natural beauty maps long strings of element labels to m_structure_dict's pointer
    std::map<std::string, std::shared_ptr<std::map<int, std::tuple<std::string, std::array<float, 3>, std::string>>>> m_element_hash_to_ptr;
    std::map<int, std::pair<std::array<float, 3>, std::array<std::array<float, 3>, 3>>> *m_transformations_dict_ptr;
    std::map<std::string, int> m_position_hash_to_atom_index;
    std::map<int, std::pair<std::string, std::array<float, 3>>> *m_structure_dict_ptr;
    std::map<int, std::tuple<std::string, std::array<float, 3>, std::string>> m_picked_structure_dict;
    std::map<int, std::tuple<std::string, std::array<float, 3>, std::string>> m_working_structure_dict;

    std::vector<std::string> m_found_element_hashes = {};

    std::array<float, 3> m_center_frac_coord;

    nlohmann::json m_json;


    public:
    explicit GeneratorFinder(std::map<int, std::pair<std::string, std::array<float, 3>>> *structure_dict_ptr, const std::string &path, std::map<int, std::pair<std::array<float, 3>, std::array<std::array<float, 3>, 3>>> *transformations_dict_ptr) {
        m_structure_dict_ptr = structure_dict_ptr;
        m_transformations_dict_ptr = transformations_dict_ptr;
        read_parse_json(std::string(path));
        build_element_hash_dict();

        // initialize m_structure_dict with the first configuration so that the cell_viewer can initialize properly
        *m_structure_dict_ptr = pop_last_column(m_element_hash_to_ptr[m_element_hash_to_ptr.begin()->first].get());
        if (m_structure_dict_ptr->size() == 0) {
            throw(std::range_error("CellViewer cannot initialize properly as m_structure_dict has size 0. Error thrown in constructor of GeneratorFinder."));
        }
    }

    void read_parse_json(std::string path) {
        // read the configuration json file whose path is "path_to_json"
        spd::info("Reading json file for parsing located at {}", path);
        std::ifstream json_file(path);
        m_json = nlohmann::json::parse(json_file);
        spd::info("Successfully parsed json");
    }

    void build_position_hash_to_atom_index() {
        // map the position hash from m_picked_structure_dict to the atom index
        for (auto &[key, val] : m_picked_structure_dict) {
            const auto &[element, frac_coord, frac_coord_hash] = val;
            m_position_hash_to_atom_index[frac_coord_hash] = key;
        }
    }

    float wrap_back_into_cell(float x) {
        return x - std::floor(x);
    }

    void apply_transformation(const std::array<float, 3> *trans_ptr, const std::array<std::array<float, 3>, 3> *rot_ptr, std::array<float, 3> *frac_coord_ptr) {
        std::array<float, 3> outer_temp{};
        for (int i = 0; i < 3; i++) {
            float inner_temp = 0;
            for (int j = 0; j < 3; j++) {
                inner_temp += (*rot_ptr)[i][j] * (*frac_coord_ptr)[j];// - m_center_frac_coord[j]);
            }
            outer_temp[i] = wrap_back_into_cell(inner_temp + (*trans_ptr)[i]);// + m_center_frac_coord[i]);
        }
        *frac_coord_ptr = outer_temp;
    }

    static std::string construct_element_hash(std::map<int, std::pair<std::string, std::array<float, 3>>> *structure_dict) {
        std::string element_hash;
        for (auto &[key, val] : *structure_dict) {
            auto [element, frac_coord] = val;
            element_hash += element;
        }
        return element_hash;
    }

    void apply_transformations(CellViewer *cell_viewer_ptr) {
        spd::info("Applying {} transformations", m_transformations_dict_ptr->size());
        m_found_element_hashes = {};
        for (auto &[transform_key, transform_val] : *m_transformations_dict_ptr) {
            if constexpr (enable_debug_logging) {
                spd::info("Showing picked transformation");
            }
            *m_structure_dict_ptr = pop_last_column(&m_picked_structure_dict);
            #if GRAPHICS
            cell_viewer_ptr->show_property_stats();
            cell_viewer_ptr->update_color_buffer();
            cell_viewer_ptr->viewer.update();
            //std::this_thread::sleep_for(std::chrono::milliseconds(500));
            #endif
            std::vector<std::string> new_position_hashes = {};
            m_working_structure_dict = m_picked_structure_dict;
            //spd::info("Before");
            //output_frac_coord(&m_working_structure_dict);
            for (auto &[key, val] : m_working_structure_dict) {
                auto &[element, frac_coord, frac_coord_hash] = val;
                apply_transformation(&transform_val.first, &transform_val.second, &frac_coord);
                // directly append transform new position to hash and emplace back
                std::string position_hash = std::format("[{:.7f},{:.7f},{:.7f}]", frac_coord[0], frac_coord[1], frac_coord[2]);
                new_position_hashes.emplace_back(position_hash);

                //spd::info("Position hash: {}", position_hash);
            }
            //spd::info("After");
            //output_frac_coord(&m_working_structure_dict);
            if constexpr (enable_debug_logging) {
                spd::info("Transformation with key {} successful", transform_key);
            }

            /*
            *m_structure_dict_ptr = pop_last_column(&m_working_structure_dict);
            if constexpr (enable_debug_logging) {
                spd::info("Showing transformed configuration!");
            }
            cell_viewer_ptr->show_property_stats();
            cell_viewer_ptr->update_color_buffer();
            cell_viewer_ptr->viewer.update(); */

            int idx = 0;
            for (const auto& new_position_hash : new_position_hashes) {
                auto &atom_index = m_position_hash_to_atom_index.at(new_position_hash);
                //spd::info("Atom index: {}", atom_index);
                auto &[element, frac_coord, frac_coord_hash]= m_picked_structure_dict.at(atom_index);
                (*m_structure_dict_ptr)[idx].first = element;
                //spd::info("Element: {}", element);
                ++idx;
            }

            // only emplace back unique element hashes
            std::string constr_element_hash = construct_element_hash(m_structure_dict_ptr);
            //spd::info("Constructed element hash: {}", constr_element_hash);
            if (m_found_element_hashes.end() == std::find(m_found_element_hashes.begin(), m_found_element_hashes.end(), constr_element_hash)) {
            //if (true) {
                m_found_element_hashes.emplace_back(constr_element_hash);
            }
            //std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }

    static std::map<int, std::pair<std::string, std::array<float, 3>>> pop_last_column(std::map<int, std::tuple<std::string, std::array<float, 3>, std::string>> *structure_dict) {
        std::map<int, std::pair<std::string, std::array<float, 3>>> temp;
        for (auto &[key, val] : *structure_dict) {
            auto [element, frac_coord, frac_coord_hash] = val;
            temp[key] = std::make_pair(element, frac_coord);
        }
        return temp;
    }

    void calculate_set_difference() {
        for (auto found_element_hash : m_found_element_hashes) {
            //spd::info("Element hash: {}", found_element_hash);
            m_element_hash_to_ptr.erase(found_element_hash);
        }
        spd::info("Found element hashes: {}", m_found_element_hashes.size());
        spd::info("Remaining element_hashes: {}", m_element_hash_to_ptr.size());
    }


    void start_reduction(CellViewer *cell_viewer_ptr) {
        bool x = true;
        int h = 0;
        do {
            spd::info("Starting reduction");
            // outer while element_hash_to_ptr.size() > 0 --> pick config
            // pick_configuration
            m_picked_structure_dict = *m_element_hash_to_ptr[m_element_hash_to_ptr.begin()->first];
            build_position_hash_to_atom_index();
            apply_transformations(cell_viewer_ptr);
            calculate_set_difference();
            // calculate_set_difference();
            x = false;
            ++h;
        //} while (x);
        } while (!m_element_hash_to_ptr.empty());
        spd::info("Finished reduction with {} found configurations", h);
    }

    void calculate_center_frac_coord() {
        float x_temp = 0;
        float y_temp = 0;
        float z_temp = 0;
        // check if all configurations have the same number of atoms, if not throw exception TO DO!!!!!!!
        auto &first_config_structure_dict = m_element_hash_to_ptr.begin()->second;
        float size = first_config_structure_dict->size();

        spd::info("Size: {}", size);
        for (auto &[key, val] : *first_config_structure_dict) {
            auto &[element, frac_coord, frac_coord_hash] = val;
            x_temp += frac_coord[0];
            y_temp += frac_coord[1];
            z_temp += frac_coord[2];
        }
        m_center_frac_coord = {x_temp/size, y_temp/size, z_temp/size};
    }

    static void output_frac_coord(std::map<int, std::tuple<std::string, std::array<float, 3>, std::string>> *structure_dict) {
        for (auto &[key, val] : *structure_dict) {
            auto &[element, frac_coord, frac_coord_hash] = val;
            std::cout << "Atom " << key << std::endl;
            for (int i= 0; i < 3; i++) {
                std::cout << frac_coord[i] << " ";
            }
            std::cout << "................" << std::endl;
        }
    }

    void build_element_hash_dict() {
        int configuration_index = 0;
        for (auto &configuration : m_json) {
            int atom_index = 0;
            std::string element_hash;
            std::map<int, std::tuple<std::string, std::array<float, 3>, std::string>> structure_dict;
            for (auto &atom : configuration) {
                element_hash += atom.at("element");
                // structure dict is comprised of atom_index (key) -> element, frac_coord, frac_coord_hash
                auto element = static_cast<std::string>(atom.at("element"));
                auto frac_coord = static_cast<std::array<float, 3>>(atom.at("frac_coord"));
                //auto frac_coord_hash = static_cast<std::string>(atom.at("frac_coord"));
                //std::string frac_coord_hash = "["+std::to_string(frac_coord[0])+","+std::to_string(frac_coord[1])+","+std::to_string(frac_coord[2])+"]";
                std::string frac_coord_hash = std::format("[{:.7f},{:.7f},{:.7f}]", frac_coord[0], frac_coord[1], frac_coord[2]);
                //spd::info(frac_coord_hash);
                structure_dict[atom_index] = std::make_tuple(element, frac_coord, frac_coord_hash);
                ++atom_index;
            }
            auto structure_dict_ptr = std::make_shared<std::map<int, std::tuple<std::string, std::array<float, 3>, std::string>>>(structure_dict);
            m_element_hash_to_ptr[element_hash] = structure_dict_ptr;
            ++configuration_index;
        }
    }
};

class GeneratorFinderFramework {
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict;
    // m_transformation is a vector (member variable) that contains .first -> translation vector and .second -> rotation matrix
    std::map<int, std::pair<std::array<float, 3>, std::array<std::array<float, 3>, 3>>> m_transformations_dict;

    GeneratorFinder m_generator_finder;
    #if GRAPHICS
    CellViewer m_cell_viewer;
    #endif

    public:
    ThreadPool thread_pool{1};
    std::queue<std::future<void>> results;

    py::dict python_transformations_dict;

    explicit GeneratorFinderFramework(const py::str &path, const py::dict &transformations)
    : m_generator_finder(&m_structure_dict, std::string(path), &m_transformations_dict)
    #if GRAPHICS
    , m_cell_viewer(&m_structure_dict)
    #endif
    {
        python_transformations_dict = transformations;
        parse_transformations();
        //check_transformation_parsing();
    }

    void parse_transformations() {
        spd::info("Parsing transformations");
        for (auto &[key, val] : python_transformations_dict) {
            m_transformations_dict[key.cast<int>()] = {val["trans"].cast<std::array<float, 3>>(), val["rot"].cast<std::array<std::array<float, 3>, 3>>()};
        }
        spd::info("Successfully parsed transformations");
    }

    void check_transformation_parsing() {
        for (auto &[key, val] : m_transformations_dict) {
            std::cout << "Translation" << std::endl;
            for (int i = 0; i < val.first.size(); i++) {
                std::cout << val.first[i] << " ";
            }
            std::cout << std::endl;

            std::cout << "Rotation" << std::endl;
            for (int i = 0; i < val.second.size(); i++) {
                for (int j = 0; j < val.second[i].size(); j++) {
                    std::cout << val.second[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }
    }
    #if GRAPHICS
    void start_reduction() {
        m_cell_viewer.initialize_viewer();
        auto m_cell_viewer_ptr = &m_cell_viewer;
        results.emplace(thread_pool.AddTask([this, m_cell_viewer_ptr]() { this->m_generator_finder.start_reduction(m_cell_viewer_ptr); }));
        m_cell_viewer.open_viewer();
    }
    #else
    void start_reduction() {
        m_generator_finder.start_reduction(nullptr);
    }
    #endif
};

void set_resource_path(std::string &path) {
	
}

PYBIND11_MODULE(finding_generators, m) {
    m.doc() = "Python bindings for functions in Python";
    #if GRAPHICS
        
    #endif

    py::class_<VirusSimulatorFramework>(m, "VirusSimulator")
    .def(py::init<py::dict, py::dict>())
    .def("dope_cell", &VirusSimulatorFramework::dope_cell)
    .def("run_simulation", &VirusSimulatorFramework::run_simulation);

    py::class_<AtomPermutatorFramework>(m, "AtomPermutator")
    .def(py::init<py::dict>())
    .def("dope_cell", &AtomPermutatorFramework::dope_cell)
    .def("run_permutation", &AtomPermutatorFramework::run_permutation);

    py::class_<GeneratorFinderFramework>(m, "GeneratorFinder")
    .def(py::init<py::str, py::dict>())
    .def("start_reduction", &GeneratorFinderFramework::start_reduction);
}
