# Sphaleron Solver

One can run the code by using the usual command lines:

gcc -o sphaleron_solve sphaleron_solve.c

./sphaleron_solve > out.dat


This repository provides a standalone C program that numerically solves the coupled radial **SU(2)‑Higgs** field equations (Eqs.&nbsp;1–5 in the accompanying manuscript) for the electroweak sphaleron.  
All major steps are fully commented inside the source file.

## Files

| File | Description |
|------|-------------|
| **sphaleron_solve_2_public.c** | Complete C implementation with detailed comments. |
| **README.md** | This file &mdash; build instructions, overview, and licensing information. |

## Build

The code uses only the C standard library and `<math.h>`. Compile with, e.g.

```bash
gcc -O2 -std=c11 sphaleron_solve_2_public.c -lm -o sphaleron_solve
```

or, for Clang:

```bash
clang -O2 -std=c11 sphaleron_solve_2_public.c -lm -o sphaleron_solve
```

## Run

```bash
./sphaleron_solve > output.dat
```

The program prints a header followed by space‑separated columns:

```
r  fA  fB  fC  H  KK  fA'  fB'  fC'  H'  KK'  vA  vB  vC  vH  vK
```

Each block of data is preceded by a comment line indicating the simulation time.

## Author

**Sarunas Verner**  
_Copyright &copy; 2025_

## License

Released under the MIT License. See the license text embedded at the top of the source file.


## Parameters

| Symbol | Default | Meaning |
|--------|---------|---------|
| `MW`   | 1.0     | \(m_W\) in units where the Higgs VEV \(v=1\) |
| `MH`   | 1.556   | \(m_H\) in the same units |

The solver automatically uses `MW^2` and `MH^2` in the radial equations, so you
can experiment with different mass ratios by simply editing these macros at the
top of **sphaleron_solve.c**.
