void run(const int cen = 0, const int seed = 32553)
{
  gSystem->Load("D0_C.so");
  D0(10000000, cen, seed);
}

