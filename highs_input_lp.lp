\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
Subject To
 R0: x0 + x1 + x2 + x3 + x6 >= 1
 R1: x0 + x1 + x2 + x4 + x7 >= 1
 R2: x0 + x1 + x2 + x5 + x8 >= 1
 R3: x0 + x3 + x4 + x5 + x6 >= 1
 R4: x1 + x3 + x4 + x5 + x7 >= 1
 R5: x2 + x3 + x4 + x5 + x8 >= 1
 R6: x0 + x3 + x6 + x7 + x8 >= 1
 R7: x1 + x4 + x6 + x7 + x8 >= 1
 R8: x2 + x5 + x6 + x7 + x8 >= 1
 R9: 4 x0 + 4 x1 + 4 x2 + 4 x3 + 4 x4 + 4 x5 + 4 x6 + 4 x7 + 4 x8 >= 8
Bounds
 x0 <= 1
 x1 <= 1
 x2 <= 1
 x3 <= 1
 x4 <= 1
 x5 <= 1
 x6 <= 1
 x7 <= 1
 x8 <= 1
End
