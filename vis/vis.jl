using JuMP
using Katana
using Plots

function visualiseKatanaCuts(m :: JuMP.Model;
                             xlims = (-2,2),
                             ylims = xlims,
                             stepsz = 0.1,
                             gif = false,
                             gifname = "anim")
    rd = 4

    katana_model = getKatanaModel(m)
    sol = MathProgBase.getsolution(katana_model)
    table = getKatanaCuts(katana_model)
    psols = getKatanaSols(katana_model) # get each solve of the LP

    xs = sol[1]
    ys = sol[2]

    # TODO allow caller to specify options, for now just plot from -5,5
    x = collect(xlims[1]:stepsz:xlims[2])

    plotlyjs()
    plt = plot(size=(800,600)) # empty plot

    M,N = size(table)
    @assert N == 5 # only concerning ourselves with 2d case rn
    anim = gif ? Animation() : nothing
    for i in 1:M
        if table[i,3] == 0 # ignore any cuts on the objective for now
            # we have a_1*x_1 + a_2*x_2 = c essentially
            # can restate x_2 as affine transform of x_1: x_2 = -(a_1/a_2)x_1 + c/a_2
            a_1 = table[i,1]
            a_2 = table[i,2]
            c = table[i,4]
            cmp = table[i,5] > 0 ? "≥" : "≤"
            yi = -(a_1/a_2)*x + c/a_2
            constr_lab = "$(round(a_1,rd))x + $(round(a_2,rd))y $cmp $(round(c,rd))"
            plot!(plt, x, yi, xlims=xlims, ylims=ylims, lab=constr_lab)
            gif && frame(anim, plt)
        end
    end

    # add scatter points of LP solutions along the way
    for xstar in psols
       scatter!(plt, [xstar[1]], [xstar[2]], lab="") 
    end

    scatter!(plt, [xs], [ys], lab="solution", color=:blue)

    gif && frame(anim, plt) # add extra frame at end
    gif && frame(anim, plt) # add extra frame at end
    gif && frame(anim, plt) # add extra frame at end
    gif && Plots.gif(anim, "$(gifname).gif", fps=1)
    gui(plt) # display plot

end
