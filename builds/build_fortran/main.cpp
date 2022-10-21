extern "C" {
    void testingfor(double* x);
}
 
int main(int argc, char* argv[]) {
    double y = 1.0;
 
    testingfor(&y);
}
