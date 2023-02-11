// clang-format off
#include <iostream>
#include <opencv2/opencv.hpp>
#include "rasterizer.hpp"
#include "global.hpp"
#include "Triangle.hpp"

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    // Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.

    // float a = rotation_angle / 180 * MY_PI;
    // model << cos(a), -sin(a), 0, 0,
    //     sin(a), cos(a), 0, 0,
    //     0, 0, 1, 0,
    //     0, 0, 0, 1;

    // return model;
    float pi = std::acos(-1);  
    float cos = std::cos(rotation_angle/180.0 * pi);
    float sin = std::sin (rotation_angle/180.0 * pi);
    Eigen::Matrix3f rodrigues;

    Eigen::Matrix4f model;
    Eigen::Vector3f zAxis(0, 0, 1);
    Eigen::Matrix3f N, temp;
    N << 0, -1, 0, 1, 0, 0, 0, 0, 0;
    temp = cos * Eigen::Matrix3f::Identity() + (1.0 - cos) * (zAxis * zAxis.transpose()) + sin * N;
    model << temp(0,0), temp(0,1), temp(0,2), 0, temp(1,0), temp(1,1), temp(1,2), 0, temp(2,0), temp(2,1), temp(2,2), 0, 0, 0, 0, 1;
    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    // Students will implement this function

    // Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();

    // // TODO: Implement this function
    // // Create the projection matrix for the given parameters.
    // // Then return it.
    
    // float n = zNear, f = zFar;

    // // 透视投影->正交投影  挤压
    // Eigen::Matrix4f Mpersp_orhtho;
    // Mpersp_orhtho << n, 0, 0, 0,
    //     0, n, 0, 0,
    //     0, 0, n + f, -n*f,
    //     0, 0, 1, 0;

    // // 正交投影->正则立方体
    //     // 将视锥信息为r,l,t,b
    // float fovY = eye_fov / 180 * MY_PI;// 角度转弧度
    // float t = tan(fovY / 2) * (-n), b = -t;// 朝向-z方向|n|
    // float r = aspect_ratio * t, l = -r;
    //     // 转换到正则立方体
    // Eigen::Matrix4f Mortho, Mtrans, Mscale;
    // Mtrans << 1, 0, 0, -(r + l) / 2,
    //     0, 1, 0, -(t + b) / 2,
    //     0, 0, 1, -(n + f) / 2,
    //     0, 0, 0, 1;
    // Mscale << 2 / (r - l), 0, 0, 0,
    //     0, 2 / (t - b), 0, 0,
    //     0, 0, 2 / (n - f), 0,
    //     0, 0, 0, 1;
    // Mortho = Mscale * Mtrans;
    
    // // 计算得到投影矩阵
    // projection = Mortho * Mpersp_orhtho;
    // return projection;
    // TODO: Copy-paste your implementation from the previous assignment.
    //     Students will implement this function
    float y = 2 * std::tan(eye_fov / 2.0 / 180.0 * std::acos(-1)) * -zNear;
    float x = y * aspect_ratio;
    Eigen::Matrix4f projection;
    Eigen::Matrix4f ortho, per2Ortho, orthoTemp;
    per2Ortho << zNear, 0, 0, 0,
                 0, zNear, 0, 0,
                 0, 0, zNear + zFar, -zNear * zFar,
                 0, 0, 1, 0;
    ortho << 2.0 / y , 0, 0, 0,
             0, 2.0 / x, 0, 0,
             0, 0, 2.0 /(zNear - zFar), 0,
             0, 0 ,0, 1;
    orthoTemp << 1, 0, 0, 0, 
                 0, 1, 0, 0,
                 0, 0, 1, -(zFar + zNear)/ 2.0,
                 0, 0, 0, 1;
    projection = ortho * orthoTemp * per2Ortho;
    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.
    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc == 2)
    {
        command_line = true;
        filename = std::string(argv[1]);
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0,0,5};


    std::vector<Eigen::Vector3f> pos
            {
                    {2, 0, -2},
                    {0, 2, -2},
                    {-2, 0, -2},
                    {3.5, -1, -5},
                    {2.5, 1.5, -5},
                    {-1, 0.5, -5}
            };

    std::vector<Eigen::Vector3i> ind
            {
                    {0, 1, 2},
                    {3, 4, 5}
            };

    std::vector<Eigen::Vector3f> cols
            {
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0}
            };

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);
    auto col_id = r.load_colors(cols);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
        cv::imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';
    }

    return 0;
}
// clang-format on