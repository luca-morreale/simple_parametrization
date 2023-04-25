
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/IO/OBJ_reader.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include "args/args.hxx"


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

namespace PMP = CGAL::Polygon_mesh_processing;


int main(int argc, char** argv) {

    args::ArgumentParser parser("Compute all geodesic measures needed");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::Positional<std::string> mesh_arg(parser, "mesh",   "Mesh on which compute geodesic distance (.off)");
    args::Positional<std::string> query_arg(parser, "query",  "Query files");
    args::Positional<std::string> out_arg(parser, "output", "Output file");
    args::Flag normalize(parser, "normalize", "Normalize geodesic distance", {'n', "normalize"});


    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string mesh_path    = args::get(mesh_arg);
    std::string queries_path = args::get(query_arg);
    std::string output_file  = args::get(out_arg);

    // read mesh
    Triangle_mesh mesh;
    std::ifstream input(mesh_path);
    input >> mesh;
    input.close();

    // if(!CGAL::IO::read_polygon_mesh(input_mesh, mesh) || !CGAL::is_triangle_mesh(mesh))
    if(!CGAL::is_triangle_mesh(mesh))
    {
        std::cerr << "Invalid input file." << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << "NF " << num_faces(mesh) <<std::endl;
    std::cout << "NV " << num_vertices(mesh) <<std::endl;

    double norm_factor = 1.0;
    if (normalize) {
        norm_factor = std::sqrt(PMP::area(mesh));
    }

    // prepare containers for query points
    IntList face_idx_s;
    IntList face_idx_t;
    BarycentricList barycentrics_s;
    BarycentricList barycentrics_t;

    // read query points
    int fidx_s, fidx_t;
    double w0_s, w1_s, w2_s;
    double w0_t, w1_t, w2_t;
    std::ifstream evaluation_file(queries_path);
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
        distances.push_back(result.first / norm_factor);

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

