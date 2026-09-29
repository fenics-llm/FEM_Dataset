# Problem 30 Handoff

Problem 30 has been redesigned as the Schäfer--Turek/DFG 2D-3 cylinder benchmark at Re = 100 while retaining residual-based VMS with SUPG, PSPG, and grad-div stabilization.

Completed:

- rewrote `question.md` with the benchmark geometry, conditions, parameters, and requested validation quantities;
- saved the official webpage as `DFG_2D-3_benchmark.html`;
- updated `mesh.md`, `manifest.csv`, `DATASET_STRUCTURE.md`, and `.gitignore`.

Reference targets are $C_{D,\max}=2.950921575$, $C_{L,\max}\approx0.47795$, and front-minus-back $\Delta p(8)\approx0.1116$.

Next: create the body-fitted mesh and `solution.py`, run mesh/time-step checks, save canonical fields and coefficient histories, and validate against these targets. No Problem 30 solver has been implemented or run yet.

Problem 34 was completed and validated in the same session; its changes are
included in the accompanying workspace commit.
