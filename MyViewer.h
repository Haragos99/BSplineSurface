// -*- mode: c++ -*-
#pragma once

#include <string>
#include <Eigen/Eigen>
#include <QGLViewer/qglviewer.h>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <QFile>
#include <QDir>
#include <QTextStream>
#include <QGLViewer/quaternion.h>
#include "Openfiler.hpp"
#include <fstream>
#include <iostream>
#include <vector>
#include "Bspline.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Smoother/JacobiLaplaceSmootherT.hh>
#include <QtGui/QKeyEvent>
#include <QtWidgets>
#include <QGLViewer/quaternion.h>
#include <map>
#include <algorithm>

using qglviewer::Vec;


class MyViewer : public QGLViewer {
    Q_OBJECT

public:
    explicit MyViewer(QWidget* parent);
    virtual ~MyViewer();
    inline double getCutoffRatio() const;
    inline void setCutoffRatio(double ratio);
    inline double getMeanMin() const;
    inline void setMeanMin(double min);
    inline double getMeanMax() const;
    inline void setMeanMax(double max);
    inline const double* getSlicingDir() const;
    inline void setSlicingDir(double x, double y, double z);
    inline double getSlicingScaling() const;
    inline void setSlicingScaling(double scaling);
    bool openMesh(const std::string& filename, bool update_view = true);
    bool openBezier(const std::string& filename, bool update_view = true);
    bool saveBezier(const std::string& filename);



    void Epsil() {
        auto dlg = std::make_unique<QDialog>(this);
        auto* hb1 = new QHBoxLayout;
        auto* vb = new QVBoxLayout;
        auto* ok = new QPushButton(tr("Ok"));
        connect(ok, SIGNAL(pressed()), dlg.get(), SLOT(accept()));
        ok->setDefault(true);
        QLabel* text;
        auto* sb_H = new QDoubleSpinBox;
        sb_H->setDecimals(4);
        sb_H->setSingleStep(0.0001);
        sb_H->setRange(0.0001, 1);
        hb1->addWidget(sb_H);
        hb1->addWidget(ok);
        vb->addLayout(hb1);
        dlg->setWindowTitle(tr("Skalar"));
        dlg->setLayout(vb);
        if (dlg->exec() == QDialog::Accepted) {
            epsilon = sb_H->value();
            update();
        }


    }

    void fullnes() {
        auto dlg = std::make_unique<QDialog>(this);
        auto* hb1 = new QHBoxLayout;
        auto* vb = new QVBoxLayout;
        auto* ok = new QPushButton(tr("Ok"));
        connect(ok, SIGNAL(pressed()), dlg.get(), SLOT(accept()));
        ok->setDefault(true);
        QLabel* text;
        auto* sb_H = new QDoubleSpinBox;
        sb_H->setDecimals(4);
        sb_H->setSingleStep(0.0001);
        sb_H->setRange(0.0, 1.0);
        hb1->addWidget(sb_H);
        hb1->addWidget(ok);
        vb->addLayout(hb1);
        dlg->setWindowTitle(tr("fullnes"));
        dlg->setLayout(vb);
        if (dlg->exec() == QDialog::Accepted) {
            bs.f = sb_H->value();
            update();
        }


    }


    bool transparent = false;

    float& getFrameSecond() { return FrameSecond; }
    double epsilon = 0.001;
    void show() {
        show_solid = !show_solid;
        update();
    }



    bool openBS(const std::string& filename, bool update_view);


signals:
    void startComputation(QString message);
    void midComputation(int percent);
    void endComputation();
    void displayMessage(const QString& message);

protected:
    virtual void init() override;
    virtual void draw() override;
    virtual void drawWithNames() override;
    virtual void postSelection(const QPoint& p) override;
    virtual void keyPressEvent(QKeyEvent* e) override;
    virtual void mouseMoveEvent(QMouseEvent* e) override;
    virtual QString helpString() const override;
private:
    struct MyTraits : public OpenMesh::DefaultTraits {
        using Point = OpenMesh::Vec3d; // the default would be Vec3f
        using Normal = OpenMesh::Vec3d;
        VertexTraits{
          OpenMesh::Vec3d original;
          double mean;              // approximated mean curvature
          std::vector<double> weigh;
          std::vector<double> distance;
          int idx_of_closest_bone;

        };
        VertexAttributes(OpenMesh::Attributes::Normal |
            OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    };
    using MyMesh = OpenMesh::TriMesh_ArrayKernelT<MyTraits>;
    using Vector = OpenMesh::VectorT<double, 3>;

    // Mesh
    void updateMesh(bool update_mean_range = true);
    void updateVertexNormals();
#ifdef USE_JET_FITTING
    void updateWithJetFit(size_t neighbors);
#endif
    void localSystem(const Vector& normal, Vector& u, Vector& v);
    double voronoiWeight(MyMesh::HalfedgeHandle in_he);
    void updateMeanMinMax();
    void updateMeanCurvature();

    // Bezier
    static void bernsteinAll(size_t n, double u, std::vector<double>& coeff);
    void generateMesh(size_t resolution);



    int _index(size_t i, size_t j, int nv) {
        return i * (nv)+j;
    }

    //Bspline
    void generateBSMesh(size_t resolution);
    void gen()
    {
        std::vector<Vec> c;
        c.push_back(Vec(0, 0, 0));
        c.push_back(Vec(1, 0, 1));
        c.push_back(Vec(3, 0, 2));
        c.push_back(Vec(4, 0, 1));
        c.push_back(Vec(5, 0, 0));
        c.push_back(Vec(5, 1, 1));
        c.push_back(Vec(5, 2, 2));
        c.push_back(Vec(5, 3, 2));
        c.push_back(Vec(5, 5, 2));
        c.push_back(Vec(5, 6, 0));
        c.push_back(Vec(3, 6, 1));
        c.push_back(Vec(2, 6, 1));
        c.push_back(Vec(1, 6, 0));
        c.push_back(Vec(0, 6, 0));
        c.push_back(Vec(0, 5, 1));
        c.push_back(Vec(0, 3, 2));
        c.push_back(Vec(0, 2, 2));
        c.push_back(Vec(0, 1, 1));
        std::vector<double> uk{ 0,0,0,0,0.5,1,1,1,1 };
        std::vector<double> vk{ 0,0,0,0,0.3,0.8,1,1,1,1 };
        std::vector<Vec> cp;
        cp.resize(30);
        int k = 0;
        /*
        for(int i =0; i < 5;i++)
        {
            for (int j = 0; j < 6; j++)
            {
                //cp.push_back(Vec(i, j, 0));
                if (i == 0 || j == 0 || i == 5 - 1
                    || j == 6 - 1) {

                    cp[i * (6) + j] = c[k];
                    k++;
                }
            }

        }
        */
        int nu = 5;
        int nv = 6;
        for (int i = 0; i <= nu - 1; i++)
        {
            int ind = _index(i, 0, nv);
            cp[ind] = c[k];
            k++;
        }

        for (int i = 1; i <= nv - 1; i++)
        {
            int ind = _index(nu - 1, i, nv);
            cp[ind] = c[k];
            k++;
        }

        for (int i = nu - 2; i >= 0; i--)
        {
            int ind = _index(i, nv - 1, nv);
            cp[ind] = c[k];
            k++;
        }
        for (int i = nv - 2; i >= 1; i--)
        {
            int ind = _index(0, i, nv);
            cp[ind] = c[k];
            k++;
        }


        model_type = ModelType::Bspline;
        bs = BSpline(3, 3, uk, vk, cp, 6 - 1, 5 - 1);
        generateBSMesh(50);
        setupCamera();
    }
    // Visualization
    void setupCamera();
    Vec meanMapColor(double d) const;
    void drawControlNet() const;

    void drawAxes() const;
    void drawAxesWithNames() const;
    static Vec intersectLines(const Vec& ap, const Vec& ad, const Vec& bp, const Vec& bd);

    // Other
    void fairMesh();

    //////////////////////
    // Member variables //
    //////////////////////

    enum class ModelType { NONE, MESH, BEZIER_SURFACE, SKELTON, INVERZ, Bspline } model_type;
    enum class SkelltonType { MAN, WRIST, ARM, FACE } skellton_type;
    // Mesh
    MyMesh mesh;

    BSpline bs;

    double distance(Vec p, Vec p1)
    {
        double len = sqrt(pow(p.x - p1.x, 2) + pow(p.y - p1.y, 2) + pow(p.z - p1.z, 2));

        return len;
    }







    float FrameSecond = 0.0;
    QHBoxLayout* hb1 = new QHBoxLayout;
    QLabel* text_ = new QLabel;
    QVBoxLayout* vBox = new QVBoxLayout;

    //std::map<int, double> faceAreaMap;
    std::vector<std::pair<int, double>> sortedVector;
    std::vector<std::pair<int, double>> finalarea;
    std::map< MyMesh::FaceHandle, int> sortedMap;

    // Custom comparator function to sort by values (double) in ascending order
    static bool sortByValue(const std::pair<int, double>& a, const std::pair<int, double>& b) {
        return a.second < b.second;
    }


    double mean(VectorMatrix& der, OpenMesh::Vec3d on) {
        
        auto norm = Vec(on.data());
        auto u = der[1][0];
        auto v = der[0][1];
        auto uu = der[2][0];
        auto vv = der[0][2];
        auto uv = der[1][1];
        auto e = u * u;
        auto f = u * v;
        auto g = v * v;
        auto l = uu * norm;
        auto m = uv * norm;
        auto n = vv * norm;
        return (n * e - 2 * m * f + l * g) / (e * g - f * f);
    }




    bool _homework = false;
    double median_of_area = 0.0;

    void homework()
    {
        _homework = true;
        int size = 0;
        for (auto f : mesh.faces()) {
            double area = mesh.calc_sector_area(mesh.halfedge_handle(f));
            sortedVector.push_back({ f.idx(), area });
            size++;
            auto d = mesh.calc_face_normal(f);
            auto r = d | d;
            MyMesh::Normal n;
            if (n < n) {}

        }
        std::sort(sortedVector.begin(), sortedVector.end(), sortByValue);
        int partSize = sortedVector.size() / 4;

        std::vector<std::pair<int, double>>::iterator begin = sortedVector.begin();
        std::vector<std::pair<int, double>>::iterator end = sortedVector.begin() + partSize;

        std::vector<std::pair<int, double>> part1(begin, end);

        begin = end;
        end = sortedVector.begin() + 2 * partSize;

        std::vector<std::pair<int, double>> part2(begin, end);

        begin = end;
        end = sortedVector.begin() + 3 * partSize;

        std::vector<std::pair<int, double>> part3(begin, end);

        std::vector<std::pair<int, double>> part4(end, sortedVector.end());

        finalarea.reserve(part2.size() + part3.size());

        finalarea.insert(finalarea.end(), part2.begin(), part2.end());

        finalarea.insert(finalarea.end(), part3.begin(), part3.end());

        median_of_area = median(sortedVector);
        for (const auto& entry : finalarea) {
            sortedMap[mesh.face_handle(entry.first)] = entry.first;
        }
    }

    double median(const std::vector<std::pair<int, double>>& numbers)
    {
        int size = numbers.size();
        if (size % 2 == 0) {
            double middle1 = numbers[size / 2 - 1].second;
            double middle2 = numbers[size / 2].second;
            return (middle1 + middle2) / 2.0;
        }
        else {
            return numbers[size / 2].second;
        }
    }


    void move(std::vector<Vec> newp, std::vector<Vec> old);

    // Bezier
    size_t degree[2];
    std::vector<Vec> control_points;

    float currentTime() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = now.time_since_epoch();
        return std::chrono::duration_cast<std::chrono::duration<float>>(duration).count();
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="edgeHandle"></param>

   
    void setupCameraBone();
    // Visualization
    double mean_min, mean_max, cutoff_ratio;
    bool show_control_points, show_solid, show_wireframe, show_skelton;
    enum class Visualization { PLAIN, MEAN, SLICING, ISOPHOTES, WEIGH, WEIGH2 } visualization;
    GLuint isophote_texture, environment_texture, current_isophote_texture, slicing_texture;
    Vector slicing_dir;
    double slicing_scaling;
    int selected_vertex;
    struct ModificationAxes {
        bool shown;
        float size;
        int selected_axis;
        Vec position, grabbed_pos, original_pos;
    } axes;
    std::string last_filename;
    std::ofstream of;
};
// idó kivonva a másikból
#include "MyViewer.hpp"
