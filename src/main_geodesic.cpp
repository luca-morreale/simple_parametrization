
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;
typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;
typedef Traits::Barycentric_coordinates BaryCoord;

typedef boost::graph_traits<Triangle_mesh> Graph_traits;
typedef Graph_traits::vertex_descriptor    vertex_descriptor;
typedef Graph_traits::vertex_iterator      vertex_iterator;
typedef Graph_traits::face_iterator        face_iterator;

typedef std::vector<int> IntList;
typedef std::vector<double> DoubleList;
typedef std::vector<BaryCoord> BarycentricList;


int main(int argc, char** argv) {

    if (argc < 4) {
        std::cerr << "ERROR: need to specify:" << std::endl <<
                     " - an obj file for the mesh" << std::endl <<
                     " - a txt file for the queries" << std::endl <<
                     " - an output txt file" << std::endl << std::endl;
        return 1;
    }

    std::string input_mesh(argv[1]);
    std::string queries(argv[2]);
    std::string output_file(argv[3]);

    // read mesh
    Triangle_mesh mesh;
    if(!CGAL::IO::read_polygon_mesh(input_mesh, mesh) || !CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "NF " << num_faces(mesh) <<std::endl;
    std::cout << "NV " << num_vertices(mesh) <<std::endl;

    // prepare containers for query points
    IntList face_idx_s;
    IntList face_idx_t;
    BarycentricList barycentrics_s;
    BarycentricList barycentrics_t;

    // read query points
    int fidx_s, fidx_t;
    double w0_s, w1_s, w2_s;
    double w0_t, w1_t, w2_t;
    std::ifstream evaluation_file(queries);
    while (evaluation_file >> fidx_s >> w0_s >> w1_s >> w2_s >> fidx_t >> w0_t >> w1_t >> w2_t) {
        face_idx_s.push_back(fidx_s);
        face_idx_t.push_back(fidx_t);
        barycentrics_s.push_back(BaryCoord({w0_s, w1_s, w2_s}));
        barycentrics_t.push_back(BaryCoord({w0_t, w1_t, w2_t}));

        if (w0_s + w1_s + w2_s - 1.0 > 1.0e-4) {
            std::cerr << "Error! One of the barycentric coordinates does not sum up to 1. \nTerminating!" << std::endl;
            return 1;
        }
        if (w0_t + w1_t + w2_t - 1.0 > 1.0e-4) {
            std::cerr << "Error! One of the barycentric coordinates does not sum up to 1. \nTerminating!" << std::endl;
            return 1;
        }
    }
    evaluation_file.close();
    std::cout << "Found " << face_idx_s.size() << " num of query points" << std::endl;

    // construct a shortest path query object and add a source point
    Surface_mesh_shortest_path shortest_paths(mesh);

    // iterate over all points
    DoubleList distances;
    for (uint i = 0; i < face_idx_s.size(); i++) {
        int fidx_t       = face_idx_t[i];
        int fidx_s       = face_idx_s[i];
        BaryCoord bary_s = barycentrics_s[i];
        BaryCoord bary_t = barycentrics_t[i];

        // create an iterator of faces
        face_iterator face_it_s = faces(mesh).first;
        // move it to the target face
        std::advance(face_it_s, fidx_s);
        // add the query point as source
        shortest_paths.add_source_point(*face_it_s, bary_s);

        // create an iterator of faces
        face_iterator face_it_t = faces(mesh).first;
        // move it to the target face
        // std::advance(face_it_t, 0);
        std::advance(face_it_t, fidx_t);
        // compute shortest path distance
        auto result = shortest_paths.shortest_distance_to_source_points(*face_it_t, bary_t);
        distances.push_back(result.first);

        // // remove query point
        shortest_paths.remove_all_source_points();
    }

    // save distance to file
    std::ofstream outstream(output_file);
    for (uint i = 0; i < distances.size() ; i++) {
        outstream << distances[i] << std::endl;
    }
    outstream.close();

    return 0;
}

