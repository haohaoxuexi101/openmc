#ifndef OPENMC_RANDOM_RAY_EXPONENTIAL_H
#define OPENMC_RANDOM_RAY_EXPONENTIAL_H

namespace openmc {

#if defined(__CUDACC__)
#define OPENMC_HOST_DEVICE __host__ __device__ __forceinline__
#else
#define OPENMC_HOST_DEVICE inline
#endif

OPENMC_HOST_DEVICE float cjosey_exponential(float tau)
{
  constexpr float c1n = -1.0000013559236386308f;
  constexpr float c2n = 0.23151368626911062025f;
  constexpr float c3n = -0.061481916409314966140f;
  constexpr float c4n = 0.0098619906458127653020f;
  constexpr float c5n = -0.0012629460503540849940f;
  constexpr float c6n = 0.00010360973791574984608f;
  constexpr float c7n = -0.000013276571933735820960f;

  constexpr float c0d = 1.0f;
  constexpr float c1d = -0.73151337729389001396f;
  constexpr float c2d = 0.26058381273536471371f;
  constexpr float c3d = -0.059892419041316836940f;
  constexpr float c4d = 0.0099070188241094279067f;
  constexpr float c5d = -0.0012623388962473160860f;
  constexpr float c6d = 0.00010361277635498731388f;
  constexpr float c7d = -0.000013276569500666698498f;

  float x = -tau;

  float den = c7d;
  den = den * x + c6d;
  den = den * x + c5d;
  den = den * x + c4d;
  den = den * x + c3d;
  den = den * x + c2d;
  den = den * x + c1d;
  den = den * x + c0d;

  float num = c7n;
  num = num * x + c6n;
  num = num * x + c5n;
  num = num * x + c4n;
  num = num * x + c3n;
  num = num * x + c2n;
  num = num * x + c1n;
  num = num * x;

  return num / den;
}

OPENMC_HOST_DEVICE float exponentialG(float tau)
{
  constexpr float d0n = 0.5f;
  constexpr float d1n = 0.176558112351595f;
  constexpr float d2n = 0.04041584305811143f;
  constexpr float d3n = 0.006178333902037397f;
  constexpr float d4n = 0.0006429894635552992f;
  constexpr float d5n = 0.00006064409107557148f;

  constexpr float d0d = 1.0f;
  constexpr float d1d = 0.6864462055546078f;
  constexpr float d2d = 0.2263358514260129f;
  constexpr float d3d = 0.04721469893686252f;
  constexpr float d4d = 0.006883236664917246f;
  constexpr float d5d = 0.0007036272419147752f;
  constexpr float d6d = 0.00006064409107557148f;

  float x = tau;

  float num = d5n;
  num = num * x + d4n;
  num = num * x + d3n;
  num = num * x + d2n;
  num = num * x + d1n;
  num = num * x + d0n;

  float den = d6d;
  den = den * x + d5d;
  den = den * x + d4d;
  den = den * x + d3d;
  den = den * x + d2d;
  den = den * x + d1d;
  den = den * x + d0d;

  return num / den;
}

OPENMC_HOST_DEVICE float exponentialG2(float tau)
{
  constexpr float g1n = -0.08335775885589858f;
  constexpr float g2n = -0.003603942303847604f;
  constexpr float g3n = 0.0037673183263550827f;
  constexpr float g4n = 0.00001124183494990467f;
  constexpr float g5n = 0.00016837426505799449f;

  constexpr float g1d = 0.7454048371823628f;
  constexpr float g2d = 0.23794300531408347f;
  constexpr float g3d = 0.05367250964303789f;
  constexpr float g4d = 0.006125197988351906f;
  constexpr float g5d = 0.0010102514456857377f;

  float x = tau;

  float num = g5n;
  num = num * x + g4n;
  num = num * x + g3n;
  num = num * x + g2n;
  num = num * x + g1n;
  num = num * x;

  float den = g5d;
  den = den * x + g4d;
  den = den * x + g3d;
  den = den * x + g2d;
  den = den * x + g1d;
  den = den * x + 1.0f;

  return num / den;
}

#undef OPENMC_HOST_DEVICE

} // namespace openmc

#endif // OPENMC_RANDOM_RAY_EXPONENTIAL_H
