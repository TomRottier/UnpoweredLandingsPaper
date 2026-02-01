export cs1, cs2, phasespace_region_cs

const inch = 96
const pt = 4 / 3
const cm = inch / 2.54
const full_page_fig_width = 8.27inch - 2 * 0.5inch # a4 with 0.5in margins
const full_page_fig_height = 11.69inch - 2 * 0.5inch # a4 with 0.5in margins
const cs1 = Makie.wong_colors()
const cs2 = colorschemes[:sandyterrain]
const cs3 = [colorschemes[:Ingres][4], colorschemes[:Hokusai1][5], colorschemes[:Tsimshian][4]]
const cs4 = colorschemes[:Hokusai1][[end, end - 3, end - 6]]

const bird = BezierPath(
    "m 27.458914,11.438595 c -1.16007,0.223045 -3.373286,1.131859 -3.755187,-0.639748 0.09679,-0.716794 0.282734,-1.1770978 0.635192,-1.7288018 2.59851,-4.528114 8.157236,-6.177948 13.071726,-5.744351 8.271465,0.713721 16.74589,2.705799 23.660965,7.5406118 C 56.156336,9.2390842 51.040117,8.1477292 45.863008,7.9353622 39.526839,7.5611322 33.407397,9.4960732 27.458914,11.438595 Z M 64.500513,4.1774562 c 4.77216,1.2593 9.814344,1.6803 14.254074,4.00004 4.14064,1.19057 8.98545,0.5019 12.46114,3.4963598 1.53753,1.01278 6.3061,4.15097 1.5246,3.90128 -3.68723,-0.67167 -7.46501,-0.29067 -11.14905,-1.05298 -4.00355,-1.11685 -8.07558,-0.94449 -12.089781,-0.0865 -16.205818,1.34863 -32.543732,4.64851 -48.823008,2.28556 -3.637257,-0.67309 -7.602254,-1.165 -10.401025,-3.83043 -3.1277703,-0.85692 -6.1432933,-1.4375 -9.2433823,-2.18083 -2.76423,-2.8957898 3.706655,-3.1480298 4.299023,-5.6325598 2.742826,-1.98247 6.2744223,-1.43654 9.1897173,-3.07927 5.214575,-1.39331005 10.735772,-1.27955005 16.104059,-1.54960005 11.434979,-0.56163 22.723603,1.43567005 33.873633,3.72898005 z",
    fit=true, flipy=true, flipx=true, bbox=Rect(-0.65, -0.5, 1, 1)
)

const arrow_force = BezierPath("m 0 0 l -5 -5 l 5 17.5 l 5 -17.5 z", fit=true)
const arrow_vel = BezierPath("m 0 1 l -3 -3 v -2 l 3 3 l 3 -3 v 2 z", fit=true)
const arrow_axis = BezierPath("m 0 6 c -1 -3 -3 -7 -5 -9 c 0 0 4 -0 5 2 c 1 -2 5 -2 5 -2 c -2 2 -4 6 -5 9 z", fit=true)
const arrow_force_polygon = Polygon([c.p for c in arrow_force.commands if !(c isa ClosePath)])
const arrow_vel_polygon = Polygon([c.p for c in arrow_vel.commands if !(c isa ClosePath)])

#= 
in arrow plot arrowhead tip at x = 1.0 and arrow shaft ends at x = 0.0
for arrowheads which extend < 0.0 either get a gap between arrowhead base and arrow shaft or arrowhead tip does not align with the specified end point of arrow in the plot - opted for latter option here
=#
const arrowhead = Polygon(Point2f[(0, -0.5), (1, 0), (0, 0.5)]) # default used in Makie
const arrowhead_axis = Point2f[
    (0.0, -0.25), (1.0, 0.0), (0.0, 0.25),
]
const arrowhead_vel = Point2f[
    (-3 / 5, -3 / 5), (-1 / 5, -3 / 5), (2 / 5, 0.0), (-1 / 5, 3 / 5), (-3 / 5, 3 / 5), (0, 0),
]
const arrowhead_force = Point2f[
    (-5 / 17.5, -5 / 17.5), (1 - 5 / 17.5, 0), (-5 / 17.5, 5 / 17.5), (0, 0),
]

const lift_colour = to_colormap(:Set1)[1]
const drag_colour = to_colormap(:Set1)[2]
const weight_colour = to_color(:black) #to_colormap(:Set1)[7]
const velocity_colour = to_color((0.7, 0.7, 0.7))

const plottheme = Theme(
    fonts=(
        regular="NewComputerModernMath-Regular",
        bold="NewComputerModern10-Bold",
        italic="NewComputerModern10-Italic",
        bolditalic="NewComputerModern10-BoldItalic",
    ),
    fontsize=12pt,
    palette=(color=cs1,),
    colormap=cs2,
    Axis=(
        xgridvisible=false, ygridvisible=false,
        xautolimitmargin=(0.0, 0.0),
    ),
    Lines=(linewidth=3,),
    Scatter=(strokewidth=1,),
    Arrows2D=(minshaftlength=0, tipwidth=6, shaftwidth=2),
)

const experimental_data = let
    dataf = joinpath(dirname(@__DIR__), "data", "mean_data.csv")
    data_all = readdlm(dataf, ',')
    data_5m = data_all[data_all[:, 1].==5, :]
    data_7m = data_all[data_all[:, 1].==7, :]
    data_9m = data_all[data_all[:, 1].==9, :]
    data_12m = data_all[data_all[:, 1].==12, :]
    m = 0.738
    S = 0.2098
    b = 1.07
    g = 9.81
    ρ = 1.23
    v̂ = nominal_speed(m, S, ρ, g)
    t̂ = nominal_time(m, S, ρ, g)
    d̂ = nominal_length(m, S, ρ)

    perch_height = 1.25 / d̂
    perch_spacing = 12 / d̂
    data = data_12m
    data_ν = @. sqrt(data[:, 5]^2 + data[:, 6]^2) / v̂
    data_θ = @. atan(data[:, 6], data[:, 5])
    data_χ = data[:, 3] / d̂
    data_ζ = data[:, 4] / d̂ .+ perch_height

    (ν=data_ν, θ=data_θ, χ=data_χ, ζ=data_ζ, perch_height, perch_spacing)
end
