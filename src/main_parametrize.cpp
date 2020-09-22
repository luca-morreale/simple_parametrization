
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <cstdlib>
#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>       Kernel;
typedef Kernel::Point_2                      Point_2;
typedef Kernel::Point_3                      Point_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>  SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor  halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor    vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor      face_descriptor;
namespace SMP = CGAL::Surface_mesh_parameterization;

typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>  UV_pmap;
typedef SMP::Square_border_uniform_parameterizer_3<SurfaceMesh>  Border_parameterizer;
typedef SMP::Barycentric_mapping_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;


bool write_obj(std::ofstream &out, SurfaceMesh & sm, UV_pmap uv_map);

int main(int argc, char** argv)
{
    std::cout << argv[1] << std::endl;
    std::ifstream in((argc>1) ? argv[1] : "data/three_peaks.off");
    if(!in) {
        std::cerr << "Problem loading the input data" << std::endl;
        return EXIT_FAILURE;
    }
    SurfaceMesh sm;
    in >> sm;
    // A halfedge on the border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;
    // The 2D points of the uv parametrisation will be written into this map
    UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("v:uv").first;
    Parameterizer param = Parameterizer(); // set corner constrain here
    SMP::Error_code err = SMP::parameterize(sm, param, bhd, uv_map);
    if(err != SMP::OK) {
        std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
    }
    
    std::ofstream out("result.obj");
    write_obj(out, sm, uv_map);

    return EXIT_SUCCESS;

}

bool write_obj(std::ofstream &out, SurfaceMesh & sm, UV_pmap uv_map)  
{
    typedef boost::unordered_map<vertex_descriptor, std::size_t> Vertex_index_map;
    
    std::size_t vertices_counter = 0, faces_counter = 0;
    
    Vertex_index_map vium;
    boost::associative_property_map<Vertex_index_map> vimap(vium);
    
    for(auto vd : sm.vertices()){
        out << "v " << sm.point(vd) << std::endl;
    }
    
    SurfaceMesh::Vertex_range::iterator  vb, ve;
    for(boost::tie(vb, ve) = sm.vertices(); vb != ve; ++vb){
        auto uv = get(uv_map, *vb);
        
        out << "vt " << (uv.x()*2.0 -1.0) << " " << (uv.y()*2.0 -1.0) << std::endl;
    }
    
    boost::graph_traits<SurfaceMesh>::vertex_iterator vit, vend;
    boost::tie(vit, vend) = vertices(sm);
    while(vit!=vend)
    {
      vertex_descriptor vd = *vit++;
      put(vimap, vd, vertices_counter++);
    }

    BOOST_FOREACH(face_descriptor fd, faces(sm)){
      halfedge_descriptor hd = halfedge(fd, sm);
      out << "f";
      BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, sm)){
        out << " " << (get(vimap, vd)+1) << "/" << (get(vimap, vd)+1) << "/" << (get(vimap, vd)+1);
      }
      out << '\n';
      faces_counter++;
    }
    if(vertices_counter != sm.number_of_vertices())
      return 0;
    else if(faces_counter != sm.number_of_faces())
      return 0;
    else
      return 1;
}


