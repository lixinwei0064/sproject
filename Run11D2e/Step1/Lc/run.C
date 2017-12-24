void run(const int cen = 0, const int seed = 32553)
{
  gSystem->Load("Lc_C.so");
  Lc(10000000, cen, seed);
}

