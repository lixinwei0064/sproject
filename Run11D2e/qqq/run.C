void run(const int cen = 3, const int seed = 32553)
{
  gSystem->Load("charm_Lc_Levy_C.so");
  charm_Lc_Levy(1e9, cen, seed);
}

