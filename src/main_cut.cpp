
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/Seam_mesh.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Surface_mesh_parameterization/Square_border_parameterizer_3.h>
#include <CGAL/Surface_mesh_parameterization/Barycentric_mapping_parameterizer_3.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Inverse_index.h>
#include <iostream>
#include <fstream>

#include <CGAL/IO/OFF_reader.h>

#include <CGAL/Surface_mesh_parameterization/internal/Containers_filler.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>


typedef CGAL::Simple_cartesian<double>      Kernel;
typedef Kernel::Point_2                     Point_2;
typedef Kernel::Point_3                     Point_3;
typedef CGAL::Polyhedron_3<Kernel>          PolyMesh;

typedef boost::graph_traits<PolyMesh>::edge_descriptor SM_edge_descriptor;
typedef boost::graph_traits<PolyMesh>::halfedge_descriptor SM_halfedge_descriptor;
typedef boost::graph_traits<PolyMesh>::vertex_descriptor SM_vertex_descriptor;

typedef CGAL::Unique_hash_map<SM_halfedge_descriptor, Point_2> UV_uhm;
typedef CGAL::Unique_hash_map<SM_edge_descriptor, bool> Seam_edge_uhm;
typedef CGAL::Unique_hash_map<SM_vertex_descriptor, bool> Seam_vertex_uhm;
typedef boost::associative_property_map<UV_uhm> UV_pmap;
typedef boost::associative_property_map<Seam_edge_uhm> Seam_edge_pmap;
typedef boost::associative_property_map<Seam_vertex_uhm> Seam_vertex_pmap;


typedef CGAL::Seam_mesh<PolyMesh, Seam_edge_pmap, Seam_vertex_pmap> SeamMesh;

typedef boost::graph_traits<SeamMesh>::vertex_descriptor   vertex_descriptor;
typedef boost::graph_traits<SeamMesh>::halfedge_descriptor halfedge_descriptor;
typedef boost::graph_traits<SeamMesh>::face_descriptor     face_descriptor;

namespace SMP = CGAL::Surface_mesh_parameterization;


typedef SMP::Square_border_uniform_parameterizer_3<SeamMesh>  Border_parameterizer;
typedef SMP::Barycentric_mapping_parameterizer_3<SeamMesh, Border_parameterizer> Parameterizer;


bool write_obj(std::ofstream &out, SeamMesh &seam_mesh, UV_pmap uv_map);


int main(int argc, char** argv)
{
    std::cout << argv[1] << std::endl;
    std::ifstream in((argc>1) ? argv[1] : "data/three_peaks.off");
    if(!in) {
        std::cerr << "Problem loading the input data" << std::endl;
        return EXIT_FAILURE;
    }
    PolyMesh sm;
    in >> sm;
    
    Seam_edge_uhm seam_edge_uhm(false);
    Seam_edge_pmap seam_edge_pm(seam_edge_uhm);
    
    Seam_vertex_uhm seam_vertex_uhm(false);
    Seam_vertex_pmap seam_vertex_pm(seam_vertex_uhm);
    
    
    const char* filename = (argc>2) ? argv[2] : "lion.selection.txt";
    SeamMesh mesh(sm, seam_edge_pm, seam_vertex_pm);
    SM_halfedge_descriptor smhd = mesh.add_seams(filename);
    if(smhd == SM_halfedge_descriptor() ) {
        std::cerr << "Warning: No seams in input" << std::endl;
    }
    // The 2D points of the uv parametrisation will be written into this map
    // Note that this is a halfedge property map, and that uv values
    // are only stored for the canonical halfedges representing a vertex
    // A halfedge on the (possibly virtual) border
    halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh, CGAL::Polygon_mesh_processing::parameters::all_default()).first;
  
    UV_uhm uv_uhm;
    UV_pmap uv_pm(uv_uhm);
        
    //// A halfedge on the (possibly virtual) border
    //halfedge_descriptor bhd = CGAL::Polygon_mesh_processing::longest_border(mesh).first;
    Parameterizer param = Parameterizer();
    SMP::parameterize(mesh, param, bhd, uv_pm);


    std::ofstream out("result.off");
    SMP::IO::output_uvmap_to_off(mesh, bhd, uv_pm, out);
    out.close();
    
    std::ofstream out2("result_mesh.obj");

    
    typedef boost::unordered_map<vertex_descriptor, std::size_t> Vertex_index_map;
    Vertex_index_map vium;
    boost::associative_property_map<Vertex_index_map> vimap(vium);
    std::size_t vertices_counter = 0, faces_counter = 0;
    
    boost::property_map<SeamMesh, CGAL::vertex_point_t>::type vpm = get(CGAL::vertex_point, mesh);
    boost::graph_traits<SeamMesh>::vertex_iterator vit, vend;
    boost::tie(vit, vend) = vertices(mesh);
    
    while(vit!=vend)
    {
        vertex_descriptor vd = *vit++;
        halfedge_descriptor hd = halfedge(vd, mesh);
        
        auto pt = get(vpm, target(hd, mesh));
        auto uv = get(uv_pm, hd);
        out2 << "v "  << pt << std::endl;
        out2 << "vt " << uv << std::endl;
        
        put(vimap, vd, vertices_counter++);
    }

    BOOST_FOREACH(face_descriptor fd, faces(mesh)) {
        halfedge_descriptor hd = halfedge(fd, mesh);
        out2 << "f";
        BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(hd, mesh)){
            //auto pt = get(vpm, source(hd, mesh));
            //pt.index();
            auto idx = get(vimap, vd) + 1;
            out2 << " " << idx << "/" << idx << "/" << idx;
        }
        out2 << std::endl;
        faces_counter++;
    }

    return EXIT_SUCCESS;

}

bool write_obj(std::ofstream &out, SeamMesh &seam_mesh, UV_pmap uv_map)
{
    
}
