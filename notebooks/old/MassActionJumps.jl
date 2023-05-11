

rateidxs = [1, 2]    # i.e. [β,ν]
reactant_stoich =
[
  [1 => 1, 2 => 1],         # 1*S and 1*I
  [2 => 1]                  # 1*I
]
net_stoich =
[
  [1 => -1, 2 => 1],        # -1*S and 1*I
  [2 => -1, 3 => 1]         # -1*I and 1*R
]
mass_act_jump = MassActionJump(reactant_stoich, net_stoich; param_idxs=rateidxs)

jump_prob = JumpProblem(prob, Direct(), mass_act_jump)
sol = solve(jump_prob, SSAStepper())