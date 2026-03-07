(* ::Package:: *)

(* ============================================================
   OCS (Bateman) — FUERZA BRUTA (composiciones) — GENERALIZADO
   + FORMATO PYTHON/SciPy + EXPORTS COMPLETOS
   ------------------------------------------------------------
   Basado en tu archivo Bateman_fuerza_bruta_final.wl:
   - Mantiene la misma lógica (preSol + chi memo + DenomFactor)
   - Mantiene tu formateador (RoundToSig/pySci/FormatPython16)
   - Generaliza para calcular TODA la cadena: X1(t)..XN(t)
   - Exporta:
       (1) Results_Bateman_Brute_FullChain_16sig.txt  (Time, X1..XN)
       (2) Results_Bateman_Brute_XN_16sig.txt         (Time, Value= XN)
   ============================================================ *)

(* -------- Precision settings -------- *)
$OCSWorkingPrecision = 32;   (* precisión interna de cálculo *)
$OCSPrintSigDigits   = 16;   (* cifras significativas de reporte *)
$MaxExtraPrecision   = 200;

toWP[x_] := SetPrecision[x, $OCSWorkingPrecision];

(* ============================================================
   Helper: si por algo llega {v,v} o lista, tomamos el primero
   ============================================================ *)
ClearAll[toScalar];
toScalar[x_] := If[ListQ[x] && Length[x] >= 1, x[[1]], x];

(* ============================================================
   1) Composiciones débiles / soluciones diofánticas (fuerza bruta)
   ============================================================ *)

ClearAll[WeakCompositions, partitionsRestricted, PrecomputeCompositions];

WeakCompositions[m_Integer?NonNegative, n_Integer?Positive] :=
  If[n == 1,
    {{m}},
    Join @@ Table[
      Prepend[#, k] & /@ WeakCompositions[m - k, n - 1],
      {k, 0, m}
    ]
  ];

partitionsRestricted[m_Integer?NonNegative, n_Integer?Positive] := WeakCompositions[m, n];

PrecomputeCompositions[maxMu_Integer?NonNegative, len_Integer?NonNegative] := Module[{tab},
  tab = Table[{}, {maxMu + 1}];
  If[len == 0,
    tab[[1]] = {{}};
    Do[tab[[u + 1]] = {}, {u, 1, maxMu}],
    Do[tab[[u + 1]] = partitionsRestricted[u, len], {u, 0, maxMu}]
  ];
  tab
];

(* ============================================================
   2) Denominador por i: Π_{j≠i} 1/( (λ_j - λ_i)^(mu_j+1) )
   ============================================================ *)

ClearAll[DenomFactor];

DenomFactor[i_Integer, Mu_List, Lambd_List] := Module[
  {li = Lambd[[i]], au2, muM},
  au2 = Delete[Lambd, i];
  muM = Delete[Mu, i];
  Product[
    1/((au2[[j]] - li)^(muM[[j]] + 1)),
    {j, 1, Length[au2]}
  ]
];

(* ============================================================
   3) OPTIMIZACIÓN PRINCIPAL: OCSPrepare
      Devuelve f(t) para el "último" nuclido del sub-sistema dado por DC
   ============================================================ *)

ClearAll[OCSPrepare];

OCSPrepare[X0_, DC_List] := Module[
  {term1, au1, mu, term2, maxMu, len, preSol, n,
   denom, auxMuTab, auxLamTab, chiMemo, chi, local},

  term1 = toWP[X0]/toWP[DC[[-1]]];

  au1 = DeleteDuplicates[DC];
  mu  = (Count[DC, #] - 1) & /@ au1;

  term2 = Product[toWP[au1[[z]]]^ (mu[[z]] + 1), {z, 1, Length[au1]}];

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
    term1 * term2 * Sum[ local[i, t] * denom[[i]], {i, 1, n} ]
  ]
];

(* ============================================================
   4) Reporte: redondeo + formato Python/SciPy
   ============================================================ *)

ClearAll[RoundToSig, pySci, FormatPython16];

RoundToSig[x_?NumericQ, n_Integer] := Module[{m, e},
  If[x == 0, 0,
    {m, e} = MantissaExponent[x, 10];
    Round[m, 10^(-(n - 1))] * 10^e
  ]
];

(* pySci: SIEMPRE notación científica estilo Python: 1.xxxe+NN *)
pySci[x_List, sig_Integer : $OCSPrintSigDigits] := pySci[toScalar[x], sig];

pySci[x_?NumericQ, sig_Integer : $OCSPrintSigDigits] := Module[
  {y = N[toScalar[x], $OCSWorkingPrecision], m, e, s, ee, dec},
  If[y == 0, "0",
    {m, e} = MantissaExponent[y, 10];

    (* Asegurar mantisa en [1,10) *)
    If[Abs[m] < 1, m = 10 m; e = e - 1];

    dec = Max[0, sig - 1];
    s = ToString @ NumberForm[m, {Infinity, dec},
      ExponentFunction -> (Null &),
      NumberPoint -> ".",
      DigitBlock -> Infinity,
      NumberSeparator -> "",
      NumberPadding -> {"", ""}
    ];

    ee = If[e >= 0, "+" <> IntegerString[e], "-" <> IntegerString[-e]];
    s <> "e" <> ee
  ]
];

FormatPython16[x_] := pySci[RoundToSig[toScalar[x], $OCSPrintSigDigits], $OCSPrintSigDigits];

(* ============================================================
   5) Generalización a TODA la cadena
   ============================================================ *)

ClearAll[ratesFromHalfLives, BuildChainFunctionsBrute, ExportFullChainBrute];

ratesFromHalfLives[halfLives_List] := (toWP[Log[2]]/toWP[#]) & /@ halfLives;

(* Devuelve lista de funciones {X1(t),...,XN(t)} *)
BuildChainFunctionsBrute[HalfLives_List, InitialCondition_] := Module[
  {N = Length[HalfLives], DC, funcs},
  funcs = Table[
    DC = ratesFromHalfLives[HalfLives[[1 ;; n]]];
    OCSPrepare[InitialCondition, DC],
    {n, 1, N}
  ];
  funcs
];

(* Exporta tabla completa: Time, X1..XN *)
ExportFullChainBrute[HalfLives_List, InitialCondition_, TimeVector_List, file_String] := Module[
  {funcs, Nn, header, rows, lines},

  funcs = BuildChainFunctionsBrute[HalfLives, InitialCondition];
  Nn = Length[funcs];

  header = StringRiffle[Join[{"Time"}, Table["X" <> ToString[k], {k, 1, Nn}]], "\t"];

  rows = Table[
    Join[
      {ToString[TimeVector[[r]]]},
      (FormatPython16[funcs[[#]][TimeVector[[r]]]] & /@ Range[Nn])
    ],
    {r, 1, Length[TimeVector]}
  ];

  lines = Join[{header}, StringRiffle[#, "\t"] & /@ rows];

  Export[file, StringRiffle[lines, "\n"], "Text"];
  Print["File '" <> file <> "' written successfully."];
];

(* ============================================================
   6) Driver (se ejecuta al hacer Get["...wl"])
   ============================================================ *)

HalfLives = {2, 2, 3, 3, 3, 4, 4, 4, 4};
InitialCondition = toWP[6.023*10^23];
TimeVector = {1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100};

(* (A) Export FULL CHAIN *)
ExportFullChainBrute[
  HalfLives,
  InitialCondition,
  TimeVector,
  "Results_Bateman_Brute_FullChain_16sig.txt"
];

(* (B) Export ONLY XN (last) *)
Module[{DC, fOCS, results, tablaFinal, lines},
  DC = ratesFromHalfLives[HalfLives];
  fOCS = OCSPrepare[InitialCondition, DC];
  results = fOCS /@ TimeVector;

  Print["Resultados (Re) estilo Python/SciPy (16 sig figs):"];
  Scan[Print[FormatPython16[#]] &, results];

  tablaFinal = Transpose[{TimeVector, FormatPython16 /@ (toScalar /@ results)}];

  lines = Join[
    {"Time\tValue"},
    (ToString[#[[1]]] <> "\t" <> #[[2]] &) /@ tablaFinal
  ];

  Export["Results_Bateman_Brute_XN_16sig.txt", StringRiffle[lines, "\n"], "Text"];
  Print["File 'Results_Bateman_Brute_XN_16sig.txt' written successfully."];
];

