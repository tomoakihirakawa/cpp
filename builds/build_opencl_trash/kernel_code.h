void kernel simple_add(global const int* A,
		       global const int* B,
		       global int* C,
		       global const int* N) {
  int ID = get_global_id(0);
  int Nthreads = get_global_size(0);
  int n = N[0];

  int ratio = (n / Nthreads);  // number of elements for each thread
  int start = ratio * ID;
  int stop  = ratio * (ID + 1);

  for (int i=start; i<stop; i++)
    C[i] = A[i] + B[i];
}
