\ LP format - for model browsing. Use MPS format to capture full model detail.
Minimize
  X0 + X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8
Subject To
 R0: X0 + X1 + X2 + X3 + X6 >= 1
 R1: X0 + X1 + X2 + X4 + X7 >= 1
 R2: X0 + X1 + X2 + X5 + X8 >= 1
 R3: X0 + X3 + X4 + X5 + X6 >= 1
 R4: X1 + X3 + X4 + X5 + X7 >= 1
 R5: X2 + X3 + X4 + X5 + X8 >= 1
 R6: X0 + X3 + X6 + X7 + X8 >= 1
 R7: X1 + X4 + X6 + X7 + X8 >= 1
 R8: X2 + X5 + X6 + X7 + X8 >= 1
 R9: X0 - X1 >= 0
 R10: X1 - X2 >= 0
 R11: X2 - X3 >= 0
 R12: X3 - X4 >= 0
 R13: X4 - X5 >= 0
 R14: X5 - X6 >= 0
 R15: X6 - X7 >= 0
 R16: X7 - X8 >= 0
Bounds
 X0 <= 1
 X1 <= 1
 X2 <= 1
 X3 <= 1
 X4 <= 1
 X5 <= 1
 X6 <= 1
 X7 <= 1
 X8 <= 1
End
