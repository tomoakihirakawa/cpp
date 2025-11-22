// cl_khr_fp64拡張機能を有効にするようコンパイラに指示
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void vector_add(__global const float *a,
                         __global const float *b,
                         __global float *c)
{
    int i = get_global_id(0);
    
    // 計算結果を一時的に保持する変数
    float temp_c = 0.0f;

    // 同じ足し算を10000回繰り返す
    for (int j = 0; j < 10000; j++) {
        temp_c += a[i] + b[i];
    }

    c[i] = temp_c;
}

__kernel void vector_add_double(__global const double *a,
                         __global const double *b,
                         __global double *c)
{
    int i = get_global_id(0);

    // 計算を10000回繰り返す
    // ループを入れないと、データ転送のオーバーヘッドで
    // GPUが遅く見えるため、計算負荷をかける
    double temp_c = 0.0;
    for (int j = 0; j < 10000; j++) {
        temp_c += a[i] + b[i];
    }
    c[i] = temp_c;
}