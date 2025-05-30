
function plotbranch!(plt1, plt2, eqb1, eqb2, cbranch)

    # H_map1 = [H(i, params, 0.0) for i in eachrow(eqb1.x)]
    # H_map2 = [H(i, params, 0.0) for i in eachrow(eqb2.x)]
    
    stable_1 = [ 
                [ eqb1.stbl[i] == 1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (eqb1.x[i, end], eqb1.x[i, 2], eqb1.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                [ eqb1.stbl[i] == 1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (eqb1.x[i, end], eqb1.x[i, 3], eqb1.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
               ]

    stable_2 = [ 
                [ eqb2.stbl[i] == 1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (eqb2.x[i, end], eqb2.x[i, 2], eqb2.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                [ eqb2.stbl[i] == 1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (eqb2.x[i, end], eqb2.x[i, 3], eqb2.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
               ]

    neutral_1 = [ 
                 [ eqb1.stbl[i] == 0 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (eqb1.x[i, end], eqb1.x[i, 2], eqb1.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                 [ eqb1.stbl[i] == 0 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (eqb1.x[i, end], eqb1.x[i, 3], eqb1.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
                ]

    neutral_2 = [ 
                 [ eqb2.stbl[i] == 0 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (eqb2.x[i, end], eqb2.x[i, 2], eqb2.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                 [ eqb2.stbl[i] == 0 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (eqb2.x[i, end], eqb2.x[i, 3], eqb2.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
                ]

    unstable_1 = [ 
                  [ eqb1.stbl[i] == -1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (eqb1.x[i, end], eqb1.x[i, 2], eqb1.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                  [ eqb1.stbl[i] == -1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (eqb1.x[i, end], eqb1.x[i, 3], eqb1.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
                 ]

    unstable_2 = [ 
                  [ eqb2.stbl[i] == -1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (eqb2.x[i, end], eqb2.x[i, 2], eqb2.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                  [ eqb2.stbl[i] == -1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (eqb2.x[i, end], eqb2.x[i, 3], eqb2.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
                 ]
    
    plot!(plt1, neutral_1[1], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt1, unstable_1[1], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt1, stable_1[1], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)
    
    plot!(plt1, neutral_2[1], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt1, unstable_2[1], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt1, stable_2[1], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)

    plot!(plt2, neutral_1[2], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt2, unstable_1[2], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt2, stable_1[2], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)
    
    plot!(plt2, neutral_2[2], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt2, unstable_2[2], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt2, stable_2[2], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)

end

function scatterbifs!(plt1, plt2, eqb1, eqb2, clp, cbp, ms, msw)

    scatter!(plt1, [ (i.x[end], i.x[2], i.x[1]) for i in eqb1.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt1, [ (i.x[end], i.x[2], i.x[1]) for i in eqb2.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)

    scatter!(plt1, [ (i.x[end], i.x[2], i.x[1]) for i in eqb1.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt1, [ (i.x[end], i.x[2], i.x[1]) for i in eqb2.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)

    scatter!(plt2, [ (i.x[end], i.x[3], i.x[4]) for i in eqb1.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt2, [ (i.x[end], i.x[3], i.x[4]) for i in eqb2.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)

    scatter!(plt2, [ (i.x[end], i.x[3], i.x[4]) for i in eqb1.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt2, [ (i.x[end], i.x[3], i.x[4]) for i in eqb2.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)

end



function plotbranch!(H, params, plt1, plt2, eqb1::EquilibriumBranch, eqb2::EquilibriumBranch, cbranch)

    Hmap1 = [H(i, params, 0.0) for i in eachrow(eqb1.x)]
    Hmap2 = [H(i, params, 0.0) for i in eachrow(eqb2.x)]
    
    stable_1 = [ 
                [ eqb1.stbl[i] == 1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 2], eqb1.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                [ eqb1.stbl[i] == 1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 3], eqb1.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
               ]

    stable_2 = [ 
                [ eqb2.stbl[i] == 1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 2], eqb2.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                [ eqb2.stbl[i] == 1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 3], eqb2.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
               ]

    neutral_1 = [ 
                 [ eqb1.stbl[i] == 0 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 2], eqb1.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                 [ eqb1.stbl[i] == 0 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 3], eqb1.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
                ]

    neutral_2 = [ 
                 [ eqb2.stbl[i] == 0 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 2], eqb2.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                 [ eqb2.stbl[i] == 0 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 3], eqb2.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
                ]

    unstable_1 = [ 
                  [ eqb1.stbl[i] == -1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 2], eqb1.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                  [ eqb1.stbl[i] == -1 && norm(eqb1.x[i, :] - eqb1.x[i+1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 3], eqb1.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
                 ]

    unstable_2 = [ 
                  [ eqb2.stbl[i] == -1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 2], eqb2.x[i, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                  [ eqb2.stbl[i] == -1 && norm(eqb2.x[i, :] - eqb2.x[i+1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 3], eqb2.x[i, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
                 ]
    
    plot!(plt1, neutral_1[1], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt1, unstable_1[1], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt1, stable_1[1], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)
    
    plot!(plt1, neutral_2[1], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt1, unstable_2[1], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt1, stable_2[1], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)

    plot!(plt2, neutral_1[2], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt2, unstable_1[2], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt2, stable_1[2], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)
    
    plot!(plt2, neutral_2[2], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt2, unstable_2[2], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt2, stable_2[2], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)

end

function plotbranch!(H, params, plt1, plt2, eqb1::PeriodicBranch, eqb2::PeriodicBranch, cbranch)

    Hmap1 = [H(i, params, 0.0) for i in eachrow(eqb1.x[:, 1, :])]
    Hmap2 = [H(i, params, 0.0) for i in eachrow(eqb2.x[:, 1, :])]
    
    stable_1 = [ 
                [ eqb1.stbl[i] == 1 && norm(eqb1.x[i, 1, :] - eqb1.x[i+1, 1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 1, 2], eqb1.x[i, 1, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                [ eqb1.stbl[i] == 1 && norm(eqb1.x[i, 1, :] - eqb1.x[i+1, 1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 1, 3], eqb1.x[i, 1, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
               ]

    stable_2 = [ 
                [ eqb2.stbl[i] == 1 && norm(eqb2.x[i, 1, :] - eqb2.x[i+1, 1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 1, 2], eqb2.x[i, 1, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                [ eqb2.stbl[i] == 1 && norm(eqb2.x[i, 1, :] - eqb2.x[i+1, 1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 1, 3], eqb2.x[i, 1, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
               ]

    neutral_1 = [ 
                 [ eqb1.stbl[i] == 0 && norm(eqb1.x[i, 1, :] - eqb1.x[i+1, 1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 1, 2], eqb1.x[i, 1, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                 [ eqb1.stbl[i] == 0 && norm(eqb1.x[i, 1, :] - eqb1.x[i+1, 1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 1, 3], eqb1.x[i, 1, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
                ]

    neutral_2 = [ 
                 [ eqb2.stbl[i] == 0 && norm(eqb2.x[i, 1, :] - eqb2.x[i+1, 1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 1, 2], eqb2.x[i, 1, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                 [ eqb2.stbl[i] == 0 && norm(eqb2.x[i, 1, :] - eqb2.x[i+1, 1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 1, 3], eqb2.x[i, 1, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
                ]

    unstable_1 = [ 
                  [ eqb1.stbl[i] == -1 && norm(eqb1.x[i, 1, :] - eqb1.x[i+1, 1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 1, 2], eqb1.x[i, 1, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ],
                  [ eqb1.stbl[i] == -1 && norm(eqb1.x[i, 1, :] - eqb1.x[i+1, 1, :]) < 1.0 ? (Hmap1[i], eqb1.x[i, 1, 3], eqb1.x[i, 1, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb1.stbl)-1 ]
                 ]

    unstable_2 = [ 
                  [ eqb2.stbl[i] == -1 && norm(eqb2.x[i, 1, :] - eqb2.x[i+1, 1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 1, 2], eqb2.x[i, 1, 1]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ],
                  [ eqb2.stbl[i] == -1 && norm(eqb2.x[i, 1, :] - eqb2.x[i+1, 1, :]) < 1.0 ? (Hmap2[i], eqb2.x[i, 1, 3], eqb2.x[i, 1, 4]) : (NaN, NaN, NaN) for i in 1:length(eqb2.stbl)-1 ]
                 ]
    
    plot!(plt1, neutral_1[1], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt1, unstable_1[1], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt1, stable_1[1], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)
    
    plot!(plt1, neutral_2[1], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt1, unstable_2[1], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt1, stable_2[1], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)

    plot!(plt2, neutral_1[2], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt2, unstable_1[2], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt2, stable_1[2], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)
    
    plot!(plt2, neutral_2[2], label = "", color = cbranch , linestyle = :solid, linewidth = 0.5)
    plot!(plt2, unstable_2[2], label = "", color = cbranch , linestyle = :dash , linewidth = 1.5)
    plot!(plt2, stable_2[2], label = "", color = cbranch , linestyle = :solid , linewidth = 1.5)

end

function scatterbifs!(H, params, plt1, plt2, eqb1, eqb2, clp, cbp, ms, msw)

    scatter!(plt1, [ (H(i.x, params, 0.0), i.x[2], i.x[1]) for i in eqb1.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt1, [ (H(i.x, params, 0.0), i.x[2], i.x[1]) for i in eqb2.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)

    scatter!(plt1, [ (H(i.x, params, 0.0), i.x[2], i.x[1]) for i in eqb1.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt1, [ (H(i.x, params, 0.0), i.x[2], i.x[1]) for i in eqb2.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)

    scatter!(plt2, [ (H(i.x, params, 0.0), i.x[3], i.x[4]) for i in eqb1.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt2, [ (H(i.x, params, 0.0), i.x[3], i.x[4]) for i in eqb2.lp ], color = clp, label = "", ms = ms, markerstrokewidth=msw)

    scatter!(plt2, [ (H(i.x, params, 0.0), i.x[3], i.x[4]) for i in eqb1.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)
    scatter!(plt2, [ (H(i.x, params, 0.0), i.x[3], i.x[4]) for i in eqb2.bp ], color = cbp, label = "", ms = ms, markerstrokewidth=msw)

end