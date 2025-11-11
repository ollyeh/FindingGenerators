#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <map>
#include <array>
#include <utility>
#include <string>
#include <stdexcept>
#include <random>
#include <algorithm>
#include <easy3d/renderer/drawable_points.h>
#include <easy3d/viewer/viewer.h>

namespace py = pybind11;

class Framework {
    std::map<int, std::pair<std::string, std::array<float, 3>>> m_structure_dict;
    inline static thread_local std::default_random_engine m_generator;
    std::vector<int> m_picked_swaps;

    public:
    py::dict python_dict;

    Framework(const py::dict &structure_dict, const py::kwargs& kwargs) {
        // error handling for keyword argument in the constructor
        // move to python excepts in the future
        if (not kwargs.contains("crystal_system")) {
            throw std::invalid_argument( "Specify keyword argument 'crystal_system'!" );
        }
        if (kwargs.contains("crystal_system") and kwargs["crystal_system"].cast<std::string>() != "cubic") {
            throw std::invalid_argument( "Crystal system must be a cubic!" );
        }

        python_dict = structure_dict;
        parse_dict();
    }

    void parse_dict() {
        // look for general, crystal structure = cubic, if not: break
        // parse atoms into map key: (atom_index, element), value: array[3] of fractional coords
        //map will be used later
        for (auto item : python_dict) {
            m_structure_dict[item.first.cast<int>()] = {item.second["element"].cast<std::string>(), item.second["frac_coord"].cast<std::array<float, 3>>()};
        }
    }

    void dope_cell(int n_swaps, std::string element) {
        std::uniform_int_distribution<int> m_uniform_int_distribution(0, m_structure_dict.size()-1);

        if (n_swaps >= m_structure_dict.size()) {
            throw std::invalid_argument("n_swaps must not be >= than the number of atoms! Don't be lazy and build the pure cell yourself...");
        }
        std::cout << "Swapping atoms" << std::endl;
        for (int i = 0; m_picked_swaps.size() < n_swaps or i>100; i++) {
            int temp = m_uniform_int_distribution(m_generator);
            //std::cout << (m_picked_swaps.end() != std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) << std::endl;
            if (m_picked_swaps.end() == std::find(m_picked_swaps.begin(), m_picked_swaps.end(), temp)) {
                m_picked_swaps.push_back(temp);
            }
        }
        for (int index : m_picked_swaps) {
            m_structure_dict[index].first = element;
            std::cout << index << " ";
        }
        std::cout << "with " << element << "..." << std::endl;
    }

    void display_cell() {
        easy3d::Viewer viewer("Atom viewer");
        easy3d::vec3 pos1 = {1, 1, 1};
        easy3d::vec3 pos2 = {2, 1, 1};

        std::vector<easy3d::vec3> atomic_positions;
        auto atoms = new easy3d::PointsDrawable("atoms");

        for (int i = 0; i < m_picked_swaps.size(); i++) {
            for (int j = 0; j < m_picked_swaps.size(); j++) {
                for (int k = 0; k < m_picked_swaps.size(); k++) {
                    atomic_positions.push_back({static_cast<float>(i), static_cast<float>(j), static_cast<float>(k)});
                }
            }
        }

        atoms->update_vertex_buffer(atomic_positions);
        atoms->set_uniform_coloring(easy3d::vec4(1.0f, 0.0f, 0.0f, 1.0f));  // r, g, b, a
        atoms->set_impostor_type(easy3d::PointsDrawable::SPHERE);
        atoms->set_point_size(20);
        viewer.add_drawable(std::shared_ptr<easy3d::PointsDrawable>(atoms));

        viewer.fit_screen();
        viewer.run();


        //initialize an instance of easy3d
        //plot atoms based on positions and label them by color
        // configure camera controls
        // if possible run on separate thread -> calling other functions will update the atomic positions (swapping can be looked at)
    }

    void generate_configurations() {
        // every cell with index in m_picked_swaps is treated as an infected virus parient
        // rest is susceptible
        //
    }
};

PYBIND11_MODULE(finding_generators, m) {
    m.doc() = "Python bindings for functions in Python";
    py::class_<Framework>(m, "Framework")
    .def(py::init<py::dict, py::kwargs>())
    .def("dope_cell", &Framework::dope_cell)
    .def("display_cell", &Framework::display_cell);
}
