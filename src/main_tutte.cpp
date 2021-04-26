
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Simple_cartesian<double>       Kernel;
typedef Kernel::Point_2                      Point_2;
typedef Kernel::Point_3                      Point_3;

typedef CGAL::Surface_mesh<Kernel::Point_3>  SurfaceMesh;
typedef CGAL::Polyhedron_3<Kernel>           PolyMesh;

typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor   vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor     face_descriptor;

typedef SurfaceMesh::Property_map<vertex_descriptor, Point_2>  UV_pmap;
typedef boost::unordered_map<vertex_descriptor, std::size_t> Vertex_index_map;

namespace SMP = CGAL::Surface_mesh_parameterization;

typedef SMP::Square_border_uniform_parameterizer_3<SurfaceMesh>  Border_parameterizer;
typedef SMP::Barycentric_mapping_parameterizer_3<SurfaceMesh, Border_parameterizer> Parameterizer;

struct Compute_area:
  public std::unary_function<const PolyMesh::Facet, double>
{
  double operator()(const PolyMesh::Facet& f) const{
    return Kernel::Compute_area_3()(
      f.halfedge()->vertex()->point(),
      f.halfedge()->next()->vertex()->point(),
      f.halfedge()->opposite()->vertex()->point() );
  }
};

void write_obj(std::ofstream &out, SurfaceMesh & sm, UV_pmap uv_map);
void check_facets_area(SurfaceMesh &mesh, UV_pmap &uv_pm, halfedge_descriptor &bhd);

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "ERROR! File name missing." << std::endl;
        return 1;
    }
    
    std::string file(argv[1]);
    std::ifstream in(argv[1]);
    if(!in) {
        std::cerr << "Problem loading the input data" << std::endl;
        return EXIT_FAILURE;
    }
    
    // read mesh
    SurfaceMesh sm;
    in >> sm;
    
    // A halfedge on the border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(sm).first;
    
    // The 2D points of the uv parametrisation will be written into this map
    UV_pmap uv_map = sm.add_property_map<vertex_descriptor, Point_2>("v:uv").first;
    
    Parameterizer param;
    if (argc > 3) {
        vertex_descriptor v1(atoi(argv[2]));
        vertex_descriptor v2(atoi(argv[3]));
        vertex_descriptor v3(atoi(argv[4]));
        vertex_descriptor v4(atoi(argv[5]));
        Border_parameterizer border_param = Border_parameterizer(v1, v2, v3, v4);
        param = Parameterizer(border_param); // set corner constrain here
    } else {
        param = Parameterizer();
    }
    
    // Parametrization
    //Parameterizer param = Parameterizer(); // set corner constrain here
    SMP::Error_code err = SMP::parameterize(sm, param, bhd, uv_map);
    
    // check parametrization is OK
    if(err != SMP::OK) {
        std::cerr << "Error: " << SMP::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
    }
    
    // save file
    std::string out_file = file.substr(0, file.size()-4) + "_tutte.obj";
    std::ofstream out(out_file);
    write_obj(out, sm, uv_map);
    
    check_facets_area(sm, uv_map, bhd);

    return EXIT_SUCCESS;

}

void write_obj(std::ofstream &out, SurfaceMesh & sm, UV_pmap uv_map)  
{
    std::size_t vertices_counter = 0;
    
    Vertex_index_map vium;
    boost::associative_property_map<Vertex_index_map> vimap(vium);
    
    // vertices
    for(auto vd : sm.vertices()){
        out << "v " << sm.point(vd) << std::endl;
    }
    
    // texture coordinates
    SurfaceMesh::Vertex_range::iterator  vb, ve;
    for(boost::tie(vb, ve) = sm.vertices(); vb != ve; ++vb){
        auto uv = get(uv_map, *vb);
        out << "vt " << (uv.x()*2.0 -1.0) << " " << (uv.y()*2.0 -1.0) << std::endl;
        
        // set vertex index 
        put(vimap, *vb, vertices_counter++);
    }
    
    // faces
    BOOST_FOREACH(face_descriptor fd, faces(sm)) {
        halfedge_descriptor hd = halfedge(fd, sm);
        out << "f";
        BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, sm)) {
            std::size_t idx = get(vimap, vd) + 1;
            out << " " << idx << "/" << idx << "/" << idx;
        }
        out << std::endl;
    }

}

void check_facets_area(SurfaceMesh &mesh, UV_pmap &uv_pm, halfedge_descriptor &bhd)
{
    std::stringstream out;
    SMP::IO::output_uvmap_to_off(mesh, bhd, uv_pm, out);

    PolyMesh tmp;
    out >> tmp;
    
    Compute_area ca;
    std::size_t num_null_faces = 0;

    for (auto it = tmp.facets_begin(); it != tmp.facets_end(); it++) {
        if (ca(*it) == 0.0) {
            num_null_faces++;
        }
    }
    
    if (num_null_faces > 0) {
        std::cerr << "WARNING: " << num_null_faces << " faces have 0 area!" << std::endl;
    }

}
