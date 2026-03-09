(* ::Package:: *)

(* ============================================================
   Benchmark driver (parametric seed repetition mode):
   81 chains, normalized initial condition, direct effective lambdas
   Combinatorial vs FrobeniusSolve

   The heavy seed is built automatically from the first three nuclides
   of each chain, repeating the block {a,b,a} a user-selected number
   of times.
   ============================================================ *)

ClearAll["Global`*"];

(* ============================================================
   0) Load verified data
   IMPORTANT: change this path to your rescued file
   ============================================================ *)

(*Folder where the data is stored*)
dataFile = "C:\\Users\\AMD RYZEN 7\\Downloads\\salvado_lise.txt";

Get[dataFile];

(* Expected from salvado_lise.txt:
   - chainsNuclidesExtended
   - lambdaEff
   - chainsLambdasExtended (optional)
*)

If[!ValueQ[chainsNuclidesExtended] || !ValueQ[lambdaEff],
  Print["ERROR: chainsNuclidesExtended or lambdaEff were not loaded."];
  Abort[];
];

(* ============================================================
   1) User controls
   ============================================================ *)

(* Repeat counts of the triad {a,b,a}. Examples:
   1 -> {a,b,a}
   2 -> {a,b,a,a,b,a}
   3 -> {a,b,a,a,b,a,a,b,a}
*)
seedRepetitionCounts = {1, 2, 4,6};

$OCSWorkingPrecision = 32;
$MaxExtraPrecision   = 200;

x0Normalized = SetPrecision[1, $OCSWorkingPrecision];
tBenchmark   = SetPrecision[1, $OCSWorkingPrecision];
repetitions  = 7;

(* ============================================================
   2) Benchmark-safe nuclear data
   ============================================================ *)

lambdaEffBenchmark = Association[lambdaEff];

(* Avoid singular prefactors for effectively stable nuclides during benchmark *)
If[KeyExistsQ[lambdaEffBenchmark, "Pd-102"] && lambdaEffBenchmark["Pd-102"] == 0,
  lambdaEffBenchmark["Pd-102"] = 1.*^-12;
];

chainsLambdasBase = (lambdaEffBenchmark /@ #) & /@ chainsNuclidesExtended;

Print["--- Benchmark input loaded successfully ---"];
Print["Number of base chains: ", Length[chainsNuclidesExtended]];
Print["Base chain lengths: ", Union[Length /@ chainsNuclidesExtended]];
Print["Base lambda vector lengths: ", Union[Length /@ chainsLambdasBase]];
Print["Selected repetition counts: ", seedRepetitionCounts];

(* ============================================================
   3) Automatic construction of repeated heavy seeds
   ============================================================ *)

ClearAll[BuildRepeatedSeed, BuildRepeatedChains, BuildLambdaChainsForRepeatCount];

(* From the first triad {a,b,a}, build r repetitions:
   r=1 -> {a,b,a}
   r=2 -> {a,b,a,a,b,a}
   etc. *)
BuildRepeatedSeed[baseTriad_List, r_Integer?Positive] := Flatten@Table[baseTriad, {r}];

(* Original chains are assumed to have the form
   {5 heavy nuclides, 6 fission products}.
   We keep the same fission-product tail and replace the heavy block
   by r repetitions of the triad formed by the first three heavy nuclides. *)
BuildRepeatedChains[chains_List, r_Integer?Positive] := Module[
  {baseTriad, tail},
  Table[
    baseTriad = chains[[k, 1 ;; 3]];
    tail      = chains[[k, 6 ;; -1]];
    Join[BuildRepeatedSeed[baseTriad, r], tail],
    {k, 1, Length[chains]}
  ]
];

BuildLambdaChainsForRepeatCount[chains_List, dict_, r_Integer?Positive] :=
  (dict /@ #) & /@ BuildRepeatedChains[chains, r];

(* ============================================================
   4) Utility functions
   ============================================================ *)

toWP[x_] := SetPrecision[x, $OCSWorkingPrecision];

(* ============================================================
   5) COMBINATORIAL SOLVER (direct lambdas, not half-lives)
   ============================================================ *)

ClearAll[
  WeakCompositions, partitionsRestricted, PrecomputeCompositions,
  DenomFactor, OCSPrepareFromRates,
  BuildChainFunctionsBruteFromRates, SolveFullChainBruteFromRates
];

WeakCompositions[m_Integer?NonNegative, n_Integer?Positive] :=
  If[n == 1,
    {{m}},
    Join @@ Table[
      Prepend[#, k] & /@ WeakCompositions[m - k, n - 1],
      {k, 0, m}
    ]
  ];

partitionsRestricted[m_Integer?NonNegative, n_Integer?Positive] :=
  WeakCompositions[m, n];

PrecomputeCompositions[maxMu_Integer?NonNegative, len_Integer?NonNegative] := Module[{tab},
  tab = Table[{}, {maxMu + 1}];
  If[len == 0,
    tab[[1]] = {{}};
    Do[tab[[u + 1]] = {}, {u, 1, maxMu}],
    Do[tab[[u + 1]] = partitionsRestricted[u, len], {u, 0, maxMu}]
  ];
  tab
];

DenomFactor[i_Integer, Mu_List, Lambd_List] := Module[
  {li = Lambd[[i]], au2, muM},
  au2 = Delete[Lambd, i];
  muM = Delete[Mu, i];
  Product[
    1/((au2[[j]] - li)^(muM[[j]] + 1)),
    {j, 1, Length[au2]}
  ]
];

OCSPrepareFromRates[X0_, DC_List] := Module[
  {term1, au1, mu, term2, maxMu, len, preSol, n,
   denom, auxMuTab, auxLamTab, chiMemo, chi, local},

  term1 = toWP[X0]/toWP[DC[[-1]]];

  au1 = DeleteDuplicates[DC];
  mu  = (Count[DC, #] - 1) & /@ au1;

  term2 = Product[toWP[au1[[z]]]^(mu[[z]] + 1), {z, 1, Length[au1]}];

  maxMu  = Max[mu];
  len    = Length[au1] - 1;
  preSol = PrecomputeCompositions[maxMu, len];

  n = Length[au1];

  denom = Table[DenomFactor[i, mu, au1], {i, 1, n}];

  auxMuTab  = Table[Delete[mu, i],  {i, 1, n}];
  auxLamTab = Table[Delete[au1, i], {i, 1, n}];

  chiMemo = <||>;

  chi[i_Integer, j_Integer] := Lookup[chiMemo, {i, j},
    chiMemo[{i, j}] = Module[{sols, auxMu, auxLam, li},
      sols   = preSol[[j + 1]];
      auxMu  = auxMuTab[[i]];
      auxLam = auxLamTab[[i]];
      li     = au1[[i]];
      Total[
        Times @@ MapThread[
          (Binomial[#1 + #2, #1] * (1/(li - #3)^#2)) &,
          {auxMu, #, auxLam}
        ] & /@ sols
      ]
    ]
  ];

  local[i_Integer, t_] := Module[{mui = mu[[i]], li = au1[[i]], tt = toWP[t]},
    Exp[-toWP[li] * tt] *
      Sum[(tt^l)/Factorial[l] * chi[i, mui - l], {l, 0, mui}]
  ];

  Function[{t},
    term1 * term2 * Sum[local[i, t] * denom[[i]], {i, 1, n}]
  ]
];

BuildChainFunctionsBruteFromRates[Rates_List, InitialCondition_] := Module[
  {nmax = Length[Rates]},
  Table[
    OCSPrepareFromRates[InitialCondition, Rates[[1 ;; n]]],
    {n, 1, nmax}
  ]
];

SolveFullChainBruteFromRates[Rates_List, InitialCondition_, t_] := Module[
  {funcs},
  funcs = BuildChainFunctionsBruteFromRates[Rates, InitialCondition];
  N[(#[t] &) /@ funcs, $OCSWorkingPrecision]
];

(* ============================================================
   6) FROBENIUS SOLVER (direct lambdas, not half-lives)
   ============================================================ *)

ClearAll[
  uniqueRatesAndMu, DiophantineS,
  ComputeF, OmegaF, DenominatorProduct, Generalization,
  XnForSubchainFrobFromRates, BuildChainFunctionsFrobFromRates,
  SolveFullChainFrobFromRates
];

uniqueRatesAndMu[lambd_List] := Module[{u},
  u = DeleteDuplicates[lambd];
  {u, (Count[lambd, #] - 1) & /@ u}
];

DiophantineS[milist_List, n_Integer] := Module[{max = Max[milist]},
  Table[
    Which[
      n == 1, If[u == 0, {{}}, {}],
      n == 2, {{u}},
      True,   FrobeniusSolve[ConstantArray[1, n - 1], u]
    ],
    {u, 0, max}
  ]
];

ComputeF[k_Integer?NonNegative, lambda_, t_] := (-t)^k * Exp[-lambda t];

OmegaF[i_, mi_, k_, a_List, m_List, preSol_List] := Module[{b, muEx, sols},
  b    = Delete[a, i];
  muEx = Delete[m, i];
  sols = preSol[[mi - k + 1]];
  Total@Table[
    Product[
      Binomial[muEx[[j]] + sol[[j]], muEx[[j]]] / (a[[i]] - b[[j]])^sol[[j]],
      {j, Length[b]}
    ],
    {sol, sols}
  ]
];

DenominatorProduct[i_, a_List, m_List] := Module[
  {ai = a[[i]], aExcluded, mExcluded},
  aExcluded = Delete[a, i];
  mExcluded = Delete[m, i];
  Times @@ MapThread[(#1 - ai)^(#2 + 1) &, {aExcluded, mExcluded}]
];

Generalization[t_, aList_, mList_, preSol_] := Module[
  {n = Length[aList], tt = toWP[t]},
  Total@Table[
    Module[{ai = aList[[i]], mi = mList[[i]], denom, suma},
      denom = DenominatorProduct[i, aList, mList];
      suma = Sum[
        ((-1)^k/k!) * ComputeF[k, ai, tt] * OmegaF[i, mi, k, aList, mList, preSol],
        {k, 0, mi}
      ];
      suma/denom
    ],
    {i, 1, n}
  ]
];

XnForSubchainFrobFromRates[RatesSub_List, InitialCondition_] := Module[
  {lambdaN, LambdU, mu, prefactor, PreSol, Xn},

  lambdaN = Last[RatesSub];
  {LambdU, mu} = uniqueRatesAndMu[RatesSub];

  prefactor = (toWP[InitialCondition]/lambdaN) *
    Product[LambdU[[k]]^(mu[[k]] + 1), {k, 1, Length[LambdU]}];

  PreSol = DiophantineS[mu, Length[LambdU]];

  Xn[t_] := prefactor * Generalization[t, LambdU, mu, PreSol];
  Xn
];

BuildChainFunctionsFrobFromRates[Rates_List, InitialCondition_] := Module[
  {nmax = Length[Rates]},
  Table[
    XnForSubchainFrobFromRates[Rates[[1 ;; n]], InitialCondition],
    {n, 1, nmax}
  ]
];

SolveFullChainFrobFromRates[Rates_List, InitialCondition_, t_] := Module[
  {funcs},
  funcs = BuildChainFunctionsFrobFromRates[Rates, InitialCondition];
  N[(#[t] &) /@ funcs, $OCSWorkingPrecision]
];

(* ============================================================
   7) Batch benchmark helpers
   ============================================================ *)

ClearAll[RunBatchBenchmark, PrintBenchmarkSummary, BenchmarkScenario, ScenarioLabel];

RunBatchBenchmark[solver_, chains_List, x0_, t_, reps_Integer?Positive] := Module[
  {warmup, times, checksum},

  Print["Warming up..."];
  warmup = solver[First[chains], x0, t];
  checksum = Total[N[warmup, 20]];

  times = Table[
    First @ AbsoluteTiming[
      checksum = Total @ Table[
        Total[N[solver[chains[[k]], x0, t], 20]],
        {k, 1, Length[chains]}
      ];
    ],
    {reps}
  ];

  <|
    "Times" -> times,
    "Median" -> Median[times],
    "Mean" -> Mean[times],
    "Min" -> Min[times],
    "Max" -> Max[times],
    "Checksum" -> checksum
  |>
];

PrintBenchmarkSummary[label_, resBrute_, resFrob_, speedup_] := Module[{},
  Print[" "];
  Print["================ BENCHMARK SUMMARY (", label, ") ================"];
  Print["Combinatorial times: ", resBrute["Times"]];
  Print["Combinatorial median: ", resBrute["Median"]];
  Print["Combinatorial mean:   ", resBrute["Mean"]];
  Print["Combinatorial min:    ", resBrute["Min"]];
  Print["Combinatorial max:    ", resBrute["Max"]];
  Print["Combinatorial checksum: ", resBrute["Checksum"]];
  Print[" "];
  Print["Frobenius times: ", resFrob["Times"]];
  Print["Frobenius median: ", resFrob["Median"]];
  Print["Frobenius mean:   ", resFrob["Mean"]];
  Print["Frobenius min:    ", resFrob["Min"]];
  Print["Frobenius max:    ", resFrob["Max"]];
  Print["Frobenius checksum: ", resFrob["Checksum"]];
  Print[" "];
  Print["Speedup (Combinatorial / FrobMedian): ", speedup];
  Print["==============================================================="];
];

ScenarioLabel[r_Integer?Positive] := Module[{heavyLen},
  heavyLen = 3 r;
  "Repeated heavy block count " <> ToString[r] <>
    " (heavy length " <> ToString[heavyLen] <> " + 6)"
];

BenchmarkScenario[label_, chains_List] := Module[
  {resBrute, resFrob, speedup},

  Print["--- Starting benchmark on scenario: ", label, " ---"];
  Print["Number of chains: ", Length[chains]];
  Print["Lambda vector lengths: ", Union[Length /@ chains]];
  Print["Benchmark time t = ", tBenchmark];
  Print["Normalized initial condition X1(0) = ", x0Normalized];
  Print["Repetitions = ", repetitions];

  resBrute = RunBatchBenchmark[
    SolveFullChainBruteFromRates,
    chains,
    x0Normalized,
    tBenchmark,
    repetitions
  ];

  resFrob = RunBatchBenchmark[
    SolveFullChainFrobFromRates,
    chains,
    x0Normalized,
    tBenchmark,
    repetitions
  ];

  speedup = resBrute["Median"]/resFrob["Median"];

  PrintBenchmarkSummary[label, resBrute, resFrob, speedup];

  <|
    "Label" -> label,
    "Brute" -> resBrute,
    "Frobenius" -> resFrob,
    "Speedup" -> speedup
  |>
];

(* ============================================================
   8) Build all selected scenarios automatically
   ============================================================ *)

scenarioAssociations = Table[
  Module[{chainsNuclides, chainsLambdas, label},
    chainsNuclides = BuildRepeatedChains[chainsNuclidesExtended, r];
    chainsLambdas  = (lambdaEffBenchmark /@ #) & /@ chainsNuclides;
    label = ScenarioLabel[r];

    Print["Prepared scenario: ", label];
    Print["  Number of chains: ", Length[chainsLambdas]];
    Print["  Lambda vector lengths: ", Union[Length /@ chainsLambdas]];

    <|"RepeatCount" -> r, "Label" -> label, "Chains" -> chainsLambdas|>
  ],
  {r, seedRepetitionCounts}
];

(* ============================================================
   9) Run all selected scenarios
   ============================================================ *)

results = BenchmarkScenario[#Label, #Chains] & /@ scenarioAssociations;

(* ============================================================
   10) Export compact comparison table
   ============================================================ *)

summaryRows = Prepend[
  Flatten[
    Table[
      {
        {
          results[[k, "Label"]], "Combinatorial",
          ToString[results[[k, "Combinatorial", "Median"]]],
          ToString[results[[k, "Combinatorial", "Mean"]]],
          ToString[results[[k, "Combinatorial", "Min"]]],
          ToString[results[[k, "Combinatorial", "Max"]]],
          ToString[results[[k, "Combinatorial", "Checksum"]]],
          ""
        },
        {
          results[[k, "Label"]], "Frobenius",
          ToString[results[[k, "Frobenius", "Median"]]],
          ToString[results[[k, "Frobenius", "Mean"]]],
          ToString[results[[k, "Frobenius", "Min"]]],
          ToString[results[[k, "Frobenius", "Max"]]],
          ToString[results[[k, "Frobenius", "Checksum"]]],
          ToString[results[[k, "Speedup"]]]
        }
      },
      {k, 1, Length[results]}
    ],
    1
  ],
  {"Scenario", "Method", "Median", "Mean", "Min", "Max", "Checksum", "Speedup_Brute_over_Frob"}
];

Export[
  "Benchmark_81Chains_ParametricSeed_Summary.tsv",
  StringRiffle[StringRiffle[#, "\t"] & /@ summaryRows, "\n"],
  "Text"
];

Print["Summary file written: Benchmark_81Chains_ParametricSeed_Summary.tsv"];
Print["--- Parametric benchmark finished ---"];
