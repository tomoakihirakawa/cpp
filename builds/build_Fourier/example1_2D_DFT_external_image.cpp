/*DOC_EXTRACT 0_FourierTransform

# DFTを使った畳み込み・相互相関の計算

```shell
sh clean
cmake -DCMAKE_BUILD_TYPE=Release ../ -DSOURCE_FILE=example1_2D_DFT_external_image.cpp
make
./example1_2D_DFT_external_image
```

２つの画像の畳み込み・相互相関をDFTを使って計算してみる．

<figure>
<img src="./img_kernel_shift.png" alt="img_kernel_shift.png" width="500px">
<figcaption>kernel image, shift image注意：縦軸が反転して表示させているので</figcaption>
</figure>

畳み込み&畳み込みと相互相関&畳み込みを行って結果を比較する．
２つの違いは$j-n$か$j+n$かだけ．
$j=0$として$\sum$の中でのインデックスの取り方をイメージすると，
$f_0g_0+f_1g_{-1}+f_2g_{-2}+\dots$というように，インデックスは互いに反対方向に進んでいくのが，畳み込み．
なので，$j=2N_{max}-1$でも

* $j$方向の畳み込みは，$j\geq 0$で値をもつ
* 一方で相互相関は，$j< 0$でも値をもつ

このことは，次のように図でイメージすると明確で，$j$によらない固定された核となる行列に対して，シフトさせるもう片方の行列
をまず$j$だけ進める．次に，相互相関の場合は，２つの行列が重なり合う要素だけで内積をとる．
畳み込みの場合は，$j$を軸にして，シフト行列を反転させてから内積をとる．
明らかに，$j<0$でも相互相関は値を持ち，畳み込みの場合は$j<0$では値を持たないものの，
行列サイズよりも大きな$j$でも値を持ち得る．

<figure>
<img src="./img_conv_corr.png" alt="img_conv_corr.png" width="500px">
<figcaption>畳み込み&畳み込み（左）と相互相関&畳み込み（右）の結果</figcaption>
</figure>

<figure>
<img src="./img_process_conv_corr.png" alt="img_process_conv_corr.png" width="500px">
<figcaption>畳み込みと相互相関の処理の捉え方</figcaption>
</figure>

<figure>
<img src="./img_process_flipped_conv.png" alt="img_process_flipped_conv.png" width="500px">
<figcaption>図を$j$方向に反転させて畳み込むと結果は相互相関と一致する</figcaption>
</figure>

*/

#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>
#include "lib_Fourier.hpp"

std::vector<std::vector<double>> read2Dcsv(const std::string& filename) {
   std::ifstream ifs(filename);
   std::vector<std::vector<double>> data;
   std::string line;
   while (std::getline(ifs, line)) {
      std::istringstream iss(line);
      std::vector<double> row;
      double value;
      while (iss >> value) {
         row.push_back(value);
         if (iss.peek() == ',')
            iss.ignore();
      }
      data.push_back(row);
   }
   ifs.close();
   return data;
};

std::vector<std::vector<double>> patch(std::vector<std::vector<double>> img_kernel_padded, const std::vector<std::vector<double>>& img_kernel) {
   for (int i = 0; i < img_kernel.size(); ++i)
      for (int j = 0; j < img_kernel[0].size(); ++j)
         img_kernel_padded[i][j] = img_kernel[i][j];
   return img_kernel_padded;
};

void write2Dcsv(const std::string& filename, const std::vector<std::vector<double>>& data) {
   std::ofstream ofs(filename);
   for (const auto& row : data) {
      for (size_t j = 0; j < row.size(); ++j) {
         ofs << std::setprecision(15) << row[j];
         if (j < row.size() - 1)
            ofs << ",";
      }
      ofs << std::endl;
   }
   std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
   ofs.close();
};

void write2Dcsv(const std::string& filename, const std::vector<std::vector<std::complex<double>>>& data) {
   std::ofstream ofs(filename);
   for (const auto& row : data) {
      for (size_t j = 0; j < row.size(); ++j) {
         ofs << std::setprecision(15) << row[j].real();
         if (j < row.size() - 1)
            ofs << ",";
      }
      ofs << std::endl;
   }
   std::cout << "write2Dcsv: " << filename << ", size = {" << data.size() << ", " << data[0].size() << "}" << std::endl;
   ofs.close();
};

int main() {

   //! read ./img_kernel.csv
   std::vector<std::vector<double>> img_kernel = read2Dcsv("./IMGkernel.csv");
   for (int j = 0; j < img_kernel[0].size(); ++j)
      std::cout << img_kernel[0][j] << ", ";
   std::vector<std::vector<double>> img_shift = read2Dcsv("./IMGshift.csv");
   write2Dcsv("./img_shift_cpp.csv", img_kernel);

   std::vector<std::vector<double>> img_kernel_padded(img_shift.size() + img_kernel.size() - 1, std::vector<double>(img_shift[0].size() + img_kernel[0].size() - 1, 0));
   std::vector<std::vector<double>> img_shift_padded(img_shift.size() + img_kernel.size() - 1, std::vector<double>(img_shift[0].size() + img_kernel[0].size() - 1, 0));
   img_kernel_padded = patch(img_kernel_padded, img_kernel);
   img_shift_padded = patch(img_shift_padded, img_shift);
   write2Dcsv("./img_kernel_pad.csv", img_kernel_padded);
   write2Dcsv("./img_shift_pad.csv", img_shift_padded);
   std::cout << "img_kernel.size() = " << img_kernel.size() << std::endl;
   std::cout << "img_shift.size() = " << img_shift.size() << std::endl;

   auto tmp = DFT(img_kernel_padded) * DFT(img_shift_padded);
   std::vector<std::vector<double>> Retmp(tmp.size(), std::vector<double>(tmp[0].size(), 0));
   for (int i = 0; i < tmp.size(); ++i)
      for (int j = 0; j < tmp[0].size(); ++j)
         Retmp[i][j] = tmp[i][j].real();
   write2Dcsv("./img_ReUV.csv", Retmp);

   auto conv = InverseDFT(tmp);  //! 問題がある
   write2Dcsv("./img_conv.csv", conv);

   /* ------------------------------------------------ */

   // Fourier2D<double> convolver2Dkernel;
   // Fourier2D<double> convolver2Dshift;
   // convolver2Dkernel.reset4convolution(img_shift.size(), img_shift[0].size(), img_kernel.size(), img_kernel[0].size());
   // convolver2Dshift.reset4convolution(img_shift.size(), img_shift[0].size(), img_kernel.size(), img_kernel[0].size());
   // convolver2Dkernel.add(img_kernel);
   // convolver2Dshift.add(img_shift, {false, false});
   // write2Dcsv("./img_conv_ByConvolver2D_CONVCONV.csv", ReConvolution(convolver2Dkernel, convolver2Dshift));

   {
      Fourier2D<double> fourier2D;
      fourier2D.reset4convolution(img_shift.size(), img_shift[0].size(), img_kernel.size(), img_kernel[0].size());
      fourier2D.add(img_shift, {false, false}, img_kernel);
      write2Dcsv("./img_conv_ByConvolver2D_CONVCONV.csv", Re(InverseDFT(fourier2D)));
   }
   {
      Fourier2D<double> fourier2D(img_shift.size(), img_shift[0].size(), img_kernel.size(), img_kernel[0].size());
      fourier2D.add(img_shift, {true, false}, img_kernel);
      // write2Dcsv("./img_shift_pad2_reversed.csv", convolver2D.shiftpadded);
      write2Dcsv("./img_conv_ByConvolver2D_CONVCORR.csv", Re(InverseDFT(fourier2D)));
   }
   // {
   //    Convolver2D<double> convolver2D;
   //    convolver2D.reset(img_shift.size(), img_shift[0].size(), img_kernel.size(), img_kernel[0].size());
   //    convolver2D.addKernel(img_kernel);
   //    convolver2D.addShift(img_shift, {true, false});
   //    convolver2D.convolve();
   //    write2Dcsv("./img_shift_pad2_reversed.csv", convolver2D.shiftpadded);
   //    write2Dcsv("./img_conv_ByConvolver2D_CONVCORR.csv", convolver2D.getReConvolution());
   // }
}
