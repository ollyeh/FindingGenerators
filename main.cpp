#include <iostream>
#include <memory>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <map>
#include <set>
#include <array>
#include <utility>
#include <string>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <easy3d/renderer/drawable.h>
#include <easy3d/renderer/drawable_lines.h>
#include <easy3d/viewer/viewer.h>
#include <thread>
#include <cstdlib>

#include "easy3d/core/point_cloud.h"
#include "easy3d/core/random.h"
#include "easy3d/renderer/drawable_points.h"
#include "easy3d/renderer/renderer.h"


namespace py = pybind11;

std::map<std::string, easy3d::vec3> colors = {
    {"Al", easy3d::vec3(1.0f, 0.0f, 0.0f)},
    {"N", easy3d::vec3(0.0f, 1.0f, 0.0f)},
    {"Sc", easy3d::vec3(0.0f, 0.0f, 1.0f)}
};


class CellViewer {
    std::vector<easy3d::vec3> m_atomic_positions;
    std::map<int, std::pair<std::string, std::array<float, 3>>> *m_ptr_structure_dict;
    std::map<std::string, easy3d::vec3> m_unique_element_colors;

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

    void open_viewer() {
        viewer.resize(100, 100);
        easy3d::PointCloud m_atoms;
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
        viewer.run();
        viewer.exit();
    }

    void close_viewer() {
        viewer.exit();
    }
};

class Agent {
    int m_state = 0;

    public:
    std::vector<int> *neighbors_ptr;

    Agent(std::vector<int> *neighbors_ptr) {
        this -> neighbors_ptr = neighbors_ptr;
    }
    void reset() {
        m_state = 0;//so 0 maps to init element before doping
    }
    void infect() {
        m_state = 1;
    }
    int get_status() {
        if (m_state == 1) {
            //reset();
        }
        return m_state;
    }
};

class Grid {
    inline static thread_local std::default_random_engine m_generator;
    inline static std::uniform_int_distribution<int> m_uniform_int_distribution;
    std::vector<std::unique_ptr<Agent>> m_agent_unique_ptrs;

public:
    explicit Grid(std::vector<int> *to_infect, std::map<int, std::vector<int>> *neighbor_dict_ptr) {
        // to infect has size 0
        std::cout << "Preparation..." << std::endl;
        m_agent_unique_ptrs.resize(neighbor_dict_ptr->size());
        for (int i = 0; i < neighbor_dict_ptr->size(); i++) {
            m_agent_unique_ptrs[i] = std::make_unique<Agent>(&neighbor_dict_ptr->at(i));
            // if agent is in to_infect vector
            if (std::find(to_infect->begin(), to_infect->end(), i) != to_infect->end()) {
                m_agent_unique_ptrs[i]->infect();
                std::cout << "Infecting agent "<< i << std::endl;
            }
            else {
                m_agent_unique_ptrs[i]->reset();
            }
        }
    }

    void update_simulation() {
        std::vector<int> agents_to_infect = {};

        for (int i = 0; i < m_agent_unique_ptrs.size(); i++) {
            // check if the current agent i is susceptible to the virus
            //std::cout << agent_infected << std::endl;
            bool const agent_infected = (m_agent_unique_ptrs[i]->get_status() == 1);

            // reference to neighbors_ptr
            auto &neighbors_ptr = m_agent_unique_ptrs[i]->neighbors_ptr;
            std::uniform_int_distribution<int> m_uniform_int_distribution(0, neighbors_ptr->size()-1);

            // pick a random neighbor from neighbors_ptr's vector that is not in temp
            int random_neighbor_index = m_uniform_int_distribution(m_generator);
            //std::cout << (m_picked_swaps.end() != std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) << std::endl;

            //chose random neighbor of agent i and check if susceptible
            int random_neighbor = neighbors_ptr->at(random_neighbor_index);

            if (agent_infected and m_agent_unique_ptrs.at(random_neighbor)->get_status() == 0 and agents_to_infect.end() == std::find(agents_to_infect.begin(), agents_to_infect.end(), random_neighbor)) {
                agents_to_infect.push_back(random_neighbor);
            }
            std::cout << "Next step we'll infect: " << agents_to_infect.size() << std::endl;
        }

        // check that agents_to_infect matches *to_infect->size()
        // think of a way of preventing the annihilation of infectious states
        // with the current logic infectious states are annihilating themselves, only one will remain after many simulation steps

        // agents have to share their next step, if it is the same agent, choose again

        for (int i = 0; i < m_agent_unique_ptrs.size(); i++) {
            bool const agent_in_agents_to_infect = agents_to_infect.end() != std::find(agents_to_infect.begin(), agents_to_infect.end(), i);
            bool const agent_infected = (m_agent_unique_ptrs[i]->get_status() == 1);
            if (agent_in_agents_to_infect) {
                //std::cout << "Agent " << i << " is infected right now." << std::endl;
                m_agent_unique_ptrs[i]->infect();
            }
            else if (agent_infected) {
                m_agent_unique_ptrs[i]->reset();
                //std::cout << "Agent " << i << " is reset right now." << std::endl;
            }
        }
    }
};


class CellOperator {
    inline static thread_local std::default_random_engine m_generator;
    std::vector<int> m_picked_swaps;
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict_init;
    std::map<int, std::pair<std::string, std::array<float, 3>>> *m_ptr_structure_dict;
    std::map<int, std::vector<int>> *m_ptr_neighbor_dict;

    public:
    CellOperator(std::map<int, std::pair<std::string, std::array<float, 3>>> *structure_dict_ptr, std::map<int, std::vector<int>> *neighbor_dict_ptr) {
        this -> m_structure_dict_init = *structure_dict_ptr;
        this ->m_ptr_structure_dict = structure_dict_ptr;
        this ->m_ptr_neighbor_dict = neighbor_dict_ptr;
    }
    //highly problematic, copy the dict and dope the copy
    //we need the undoped dict for the virus simulation
    void dope_cell(int n_swaps, std::string element) {
        std::uniform_int_distribution<int> m_uniform_int_distribution(0, m_ptr_structure_dict->size()-1);
        if (n_swaps >= m_ptr_structure_dict->size()) {
            throw std::invalid_argument("n_swaps must not be >= than the number of atoms! Don't be lazy and build the pure cell yourself...");
            }
        std::cout << "Swapping atoms" << std::endl;
        for (int i = 0; m_picked_swaps.size() < n_swaps; i++) {
            int temp = m_uniform_int_distribution(m_generator);
            //std::cout << (m_picked_swaps.end() != std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) << std::endl;
            if (m_picked_swaps.end() == std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) {
                m_picked_swaps.push_back(temp);
                }
            }
        for (int index : m_picked_swaps) {
            m_ptr_structure_dict->at(index).first = element;
            std::cout << index << " ";
        }
        std::cout << "with " << element << "..." << std::endl;
    }

    void prepare_simulation() {
        // to_indext = picked_swaps
        std::cout << m_picked_swaps.size() << std::endl;
        auto grid = Grid(&m_picked_swaps, m_ptr_neighbor_dict);
        for (int i = 0; i < 200; i++) {
            std::cout << "Simulation step: " << i << std::endl;
            grid.update_simulation();
        }



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

    void close_viewer() {
        m_cell_viewer.close_viewer();
    }

    void prepare_simulation() {
        m_cell_operator.prepare_simulation();
    }
};

PYBIND11_MODULE(finding_generators, m) {
    m.doc() = "Python bindings for functions in Python";
    py::class_<Framework>(m, "Framework")
    .def(py::init<py::dict, py::dict, py::kwargs>())
    .def("dope_cell", &Framework::dope_cell)
    .def("display_cell", &Framework::display_cell)
    .def("close_viewer", &Framework::close_viewer)
    .def("prepare_simulation", &Framework::prepare_simulation);
}
