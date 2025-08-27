#include "MyViewer.h"






bool MyViewer::openBS(const std::string& filename, bool update_view)
{

    size_t du, dv;
    size_t nu, nv;
    std::vector<double> knot_u;
    std::vector<double> knot_v;

    std::ifstream f(filename.c_str());
    f.exceptions(std::ios::failbit | std::ios::badbit);
    f >> du >> dv;
    f >> nu >> nv;
    knot_u.resize(du + nu+1);
    knot_v.resize(dv + nv+1);
    for (int i = 0; i < knot_u.size(); i++)
    {
        f >> knot_u[i];
    }
    for (int i = 0; i < knot_v.size(); i++)
    {
        f >> knot_v[i];
    }

    std::vector<Vec>cp;
    cp.resize(nu * nv);
    float x, y, z;
    for (int i = 0; i <= nu - 1; i++)
    {
        f >> x >> y >> z;
        int ind = _index(i, 0,nv);
        cp[ind] = Vec(x, y, z);

    }

    for (int i = 1; i <= nv - 1; i++)
    {
        f >> x >> y >> z;
        int ind = _index(nu - 1, i,nv);
        cp[ind] = Vec(x, y, z);

    }

    for (int i = nu - 2; i >= 0; i--)
    {
        f >> x >> y >> z;
        int ind = _index(i, nv - 1,nv);
        cp[ind] = Vec(x, y, z);

    }
    for (int i = nv - 2; i >= 1; i--)
    {
        f >> x >> y >> z;
        int ind = _index(0, i,nv);
        cp[ind] = Vec(x, y, z);

    }


    bs =BSpline(du, dv, knot_u, knot_v, cp, nv - 1, nu - 1);
    //bs.open(filename);
    model_type = ModelType::Bspline;
    updateMesh(true);
    setupCamera();
    
   return true;
}


bool MyViewer::openBezier(const std::string& filename, bool update_view) {
    size_t n, m;
    try {
        std::ifstream f(filename.c_str());
        f.exceptions(std::ios::failbit | std::ios::badbit);
        f >> n >> m;
        degree[0] = n++; degree[1] = m++;
        control_points.resize(n * m);
        for (size_t i = 0, index = 0; i < n; ++i)
            for (size_t j = 0; j < m; ++j, ++index)
                f >> control_points[index][0] >> control_points[index][1] >> control_points[index][2];
    }
    catch (std::ifstream::failure&) {
        return false;
    }
    model_type = ModelType::BEZIER_SURFACE;
    last_filename = filename;
    updateMesh(update_view);
    if (update_view)
        setupCamera();
    return true;
}


bool MyViewer::openMesh(const std::string& filename, bool update_view) {
    if (!OpenMesh::IO::read_mesh(mesh, filename) || mesh.n_vertices() == 0)
        return false;
    model_type = ModelType::MESH;
    last_filename = filename;
    updateMesh(update_view);

    if (update_view)
        setupCamera();
    
    for (auto v : mesh.vertices())
    {
        mesh.data(v).original = mesh.point(v);
    }
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();

    return true;
}



