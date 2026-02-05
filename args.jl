
const defFloat = Float64
const defInt = Int64

using Plots
using DelimitedFiles
using Measures
using LaTeXStrings

using StatsBase

# Define a custom color palette and figure defaults.
colorlist = [
    :indianred,       # CD5C5C
    :cornflowerblue,  # 6699FF
    :olivedrab,       # 6B8E23
    :darkgoldenrod,
    :mediumpurple,    # 9370DB
]
default(
    guidefont=font(10, "Computer Modern"),  # axis labels
    tickfont=font(10, "Computer Modern"),   # tick labels
    legendfont=font(10, "Computer Modern"), # legend
    titlefont=font(10, "Computer Modern"),  # title
    palette=colorlist,
    # left_margin=5pt, right_margin=5pt, bottom_margin=5pt, top_margin=5pt,
)

# Simple plot-saving helper function.
function saveplot(plt, filename::String; background=:transparent, args...)
    plot!( plt; background=background, args... )
    savefig( plt,  filename )
    plot!( plt; background=:white )
    return plt
end
