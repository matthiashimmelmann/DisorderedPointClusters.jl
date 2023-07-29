using Plots, GLMakie, LinearAlgebra, LaTeXStrings
#TODO NLayers=8 in grey
fluctuations6 = [[1.03183000008,0.76080118787,1.85546258799,0.97469608812,1.5296059698,0.83367645577], [0.81212039822,0.83035670537,0.86667582948,0.65208841456,0.65354209216,1.0334402561], 
                [0.71227975548,0.74590495977,0.66315273782,0.87108144911,0.74692216259,0.60883546638], [0.84301359351,1.03674849745,0.73845685414,0.79102132094,0.74861341318,0.83754962849],
                [1.0475818174,0.82931338792,0.94715796313,0.81664802472,0.86507047759,0.76834683832], [0.8130267525,0.78354318355,0.86361511899,0.8282565849930,0.74295853401,0.89651791784],
                [0.76009350466,0.58052938,0.63448479597,0.86286803327,0.83547264102,0.86626284814], [0.78045476105,0.85095005541,0.83552983556,0.77063341018,0.78132696108,0.76600699102]]
isotropy6 = [[2.34294056024,2.32703870439,1.22762247839,2.62516312135,1.73229067735,2.26166743066], [1.87829718424,1.29026726643,1.72962543465,1.43820203122,1.61771736216,1.36732483488], 
            [1.33719270644,1.42191814893,1.33991419472,1.65513823029,1.28679935719,1.34280308775], [1.33399876287,1.67026455548,1.36813987934,1.2191040306,1.15419509496,1.27604318603], 
            [1.56707670459,1.15557842881,1.54135185953,1.54221706805,1.31678963292,1.2030669472], [1.2110267296,1.17397733004,1.34552116772,1.14679618312,1.12599837125,1.32378228479],
            [1.27739087268,1.1475590248,1.47241186829,1.23793004622,1.22440231868,1.18889003487], [1.26932253199,1.73250629934,1.2881119598,1.51757679282,1.1905399676,1.40172924077]]
isotropy6 = [1 ./ isotropy for isotropy in isotropy6]
totalFluctuations6 = vcat(fluctuations6...)
totalIsotropy6 = vcat(isotropy6...)
mean6 = [sum(totalIsotropy6)/length(totalIsotropy6), sum(totalFluctuations6)/length(totalFluctuations6)]
X6 = [[iso-mean6[1] for iso in totalIsotropy6] [fluctu-mean6[2] for fluctu in totalFluctuations6]]
ellipse_vert_X6 = [GLMakie.Point2f0(mean6+svd(X6'*X6/length(totalIsotropy6)).V*[sqrt(1.5*svd(X6'*X6/length(totalIsotropy6)).S[1])*cos(phi), sqrt(1.5*svd(X6'*X6/length(totalIsotropy6)).S[2])*sin(phi)]) for phi in 0:0.1:2pi]
PCA6_1 = svd(X6'*X6/length(totalIsotropy6)).V[1,:] * sqrt(svd(X6'*X6/length(totalIsotropy6)).S[1])
PCA6_2 = svd(X6'*X6/length(totalIsotropy6)).V[2,:] * sqrt(svd(X6'*X6/length(totalIsotropy6)).S[2])


fluctuations4 = [[0.21336048558, 0.2123044435, 0.55280451076, 0.21140577007, 0.39355073809, 0.56426831228], [0.69227535183, 0.93918801715, 0.66541699724, 0.88723638936, 0.86145388567, 0.88957038679], 
                [0.39289538933, 0.21127587012, 0.40401427984, 0.21981433419, 0.44464272984, 0.31374216504], [0.67290107274, 0.45843118039, 0.41920654549, 0.7941841463, 0.67245863767, 0.40280434322], 
                [0.57027783927, 0.75041376834, 0.62026616015, 0.75249251886, 0.61419070613, 0.66829418765], [0.66485419012, 0.74505214904, 0.86795761265, 0.70578359241, 0.52006510357, 0.56398643414],
                [0.87241494869, 0.92690032539, 0.76412013786, 0.73529866606, 0.92305849271, 0.76300210044], [0.92028493114, 0.95462008533, 0.79287404565, 0.90492071785, 0.71913051567, 0.64124653055]]
isotropy4 = [[1.00309352991, 1.00125874934, 1.17365030981, 1.00326597086, 1.28954075948, 1.18561940852], [2.06520343035, 2.18730040723, 1.22806409681, 1.67322881174, 1.81641005152, 1.69420156603], 
            [1.53980001403, 1.00654908893, 1.51588176012, 1.00496095451, 1.18472485012, 1.02801705472], [1.54287798477, 1.08858196968, 1.11487435181, 1.34217040012, 1.14759171821, 1.33057811367], 
            [1.17209340626, 1.17441133041, 1.39411309576, 1.29788701371, 1.27406998333, 1.23994636769], [1.23994636769, 1.25187075851, 1.40078763884, 1.26367240159, 1.23113688812, 1.33493830941],
            [1.33000252677, 1.46000782151, 1.39171764453, 1.17934595283, 1.4689270229, 1.38706591967], [1.1249101683, 1.39109888914, 1.17042380966, 1.26641455797, 1.24646920407, 1.26051433235]]
isotropy4 = [1 ./ isotropy for isotropy in isotropy4]
totalFluctuations4 = vcat(fluctuations4...)
totalIsotropy4 = vcat(isotropy4...)
mean4 = [sum(totalIsotropy4)/length(totalIsotropy4), sum(totalFluctuations4)/length(totalFluctuations4)]
X4 = [[iso-mean4[1] for iso in totalIsotropy4] [fluctu-mean4[2] for fluctu in totalFluctuations4]]
ellipse_vert_X4 = [GLMakie.Point2f0(mean4+svd(X4'*X4/length(totalIsotropy4)).V*[sqrt(1.5*svd(X4'*X4/length(totalIsotropy4)).S[1])*cos(phi), sqrt(1.5*svd(X4'*X4/length(totalIsotropy4)).S[2])*sin(phi)]) for phi in 0:0.1:2pi]
PCA4_1 = svd(X4'*X4/length(totalIsotropy4)).V[1,:] * sqrt(svd(X4'*X4/length(totalIsotropy4)).S[1])
PCA4_2 = svd(X4'*X4/length(totalIsotropy4)).V[2,:] * sqrt(svd(X4'*X4/length(totalIsotropy4)).S[2])

fluctuations2 = [[0.43571270647, 0.44873678269, 0.44735820062, 0.44562079156, 0.433161593905, 0.44784838698], [0.68502490225, 1.41422649301, 1.77364311495, 1.78430149495, 1.19308141361, 1.30710285451],
                [0.21776292979, 0.21386267191, 0.54886466978, 0.21401761627, 0.21397237319, 0.21285922389], [0.68559818825, 0.69492265995, 1.25521419424, 0.89013967655, 0.73780812074, 0.75827302668],
                [0.62860962957, 1.439429212, 1.25474121725, 0.79977350763, 1.3856453155, 1.08877100903], [0.79955857562, 0.60957402489, 0.75595053191, 0.80604281658, 0.95544076834, 0.87815029573],
                [0.54192480461, 0.42092561751, 0.92343961594, 0.90004901159, 1.00999302969, 0.5150763662], [0.5973486303, 0.84965507045, 0.57413036896, 0.97729873788, 0.45092228182, 0.74591358765]]
isotropy2 = [[1.78904130549, 1.78986126855, 1.78879260942, 1.78711825823, 1.78811724, 1.78824585561], [1.59942169551, 2.89047669863, 1.96042208755, 1.96636061525, 2.58419593116, 2.77956273592],
            [1.00093678934, 1.00082760906, 1.20168602712, 1.00132862285, 1.00235296248, 1.00063747091], [1.03168129366, 1.10402071678, 1.20938234025, 1.4052953058, 1.08845176207, 1.14199935132],
            [1.34294185801, 1.66716535411, 1.64778285551, 2.27096446916, 2.16478432155, 1.08308868512], [1.30863257344, 1.18401536087, 1.1791279253, 1.51368701248, 1.45663060036, 1.46388002497],
            [1.0907210088, 1.60250859193, 1.17750132175, 1.40326512765, 1.33256213355, 1.63758704979], [1.43837675775, 1.51477684061, 1.09021617224, 1.2935833899, 1.19329727414, 1.56330950118]]
isotropy2 = [1 ./ isotropy for isotropy in isotropy2]
totalFluctuations2 = vcat(fluctuations2...)
totalIsotropy2 = vcat(isotropy2...)
mean2 = [sum(totalIsotropy2)/length(totalIsotropy2), sum(totalFluctuations2)/length(totalFluctuations2)]
X2 = [[iso-mean2[1] for iso in totalIsotropy2] [fluctu-mean2[2] for fluctu in totalFluctuations2]]
PCA2_1 = svd(X2'*X2/length(totalIsotropy2)).V[1,:] * sqrt(svd(X2'*X2/length(totalIsotropy2)).S[1])
PCA2_2 = svd(X2'*X2/length(totalIsotropy2)).V[2,:] * sqrt(svd(X2'*X2/length(totalIsotropy2)).S[2])
indices = [1,2,3,4,5,6]
ellipse_vert_X2 = [GLMakie.Point2f0(mean2+svd(X2'*X2/length(totalIsotropy2)).V*[sqrt(1.5*svd(X2'*X2/length(totalIsotropy2)).S[1])*cos(phi), sqrt(1.5*svd(X2'*X2/length(totalIsotropy2)).S[2])*sin(phi)]) for phi in 0:0.1:2pi]

h_isotropy =[3.1220698877,2.88061367421,2.66130891198,2.25467367875,1.96079109972,1.67138294638,1.39795888949,1.26234888971,1.12179671452,1.00165047221,1.03261319987,1.12842130852,1.2175453934,1.35308608101,1.4609356218,1.52477762086,1.58085699917,1.68694077341]
h_fluctuation = [1.73677904191, 1.55819996932, 1.4103094609, 1.11292295061, 0.89578347854, 0.70865795955, 0.53791899168, 0.46136207657,0.39946546196,0.3560685774,0.35409733015,0.34126142946,0.33535473501,0.34650691666,0.36450256785,0.37959595259,0.39172260532,0.42778106229]
h_isotropy = [1 ./ isotropy for isotropy in h_isotropy]

#=
plt1 = plot(xlims=(0.9,6.1), ylims=(0.15,1.5), title="Curvature Fluctuations")
foreach(t->plot!(plt1, indices, fluctuations4[t], linewidth=6, label="$(t+1) Necks", color=cgrad(:turbo)[Int(round(256*t/length(fluctuations2)))]), 1:length(fluctuations2))
#scatter!(plt1, [1,2,4,5,6], [0.21, 0.21, 0.21, 0.21, 0.21], markershape=:star, markersize=16, label="D", color=:black)
scatter!(plt1, [1,2,4], [0.21, 0.21, 0.21, 0.21, 0.21], markershape=:star, markersize=16, label="D", color=:black)

plt2 = plot(xlims=(0.9,6.1), ylims=(0.45,1.05), title="Isotropy (largest/smallest)")
foreach(t->plot!(plt2, indices, isotropy4[t], linewidth=6, label="$(t+1) Necks", color=cgrad(:turbo)[Int(round(256*t/length(fluctuations2)))]), 1:length(isotropy2))
#scatter!(plt2, [1,2,4,5,6], [1, 1, 1, 1, 1], markershape=:star, markersize=16, label="D", color=:black)
scatter!(plt2, [1,2,4], [1, 1, 1, 1, 1], markershape=:star, markersize=16, label="D", color=:black)

plt = plot(plt1,plt2, layout=(1,2), size=(1300,800), xtickfontsize=16,ytickfontsize=16, legendfontsize=13)
savefig(plt, "FluctuationsPlot2Layers.png")
display(plt)=#

plt = GLMakie.Figure(resolution=(1050,1050),fontsize=24)
ax = Axis(plt[1,1], xlabel = L"$\beta_1^{0,2}$", ylabel = L"$(\Delta K / \Gamma^2)^2$")

GLMakie.mesh!(ax,ellipse_vert_X2; color=RGBA{Float64}(0, 1, 0, 0.15),strokewidth=0)
GLMakie.mesh!(ax,ellipse_vert_X4; color=RGBA{Float64}(1, 0, 0, 0.15),strokewidth=0)
GLMakie.mesh!(ax,ellipse_vert_X6; color=RGBA{Float64}(0, 0, 1, 0.15),strokewidth=0)

foreach(i->GLMakie.linesegments!(ax,[h_isotropy[i],h_isotropy[i+1]],[h_fluctuation[i],h_fluctuation[i+1]]; linewidth=5, color=:black), 1:length(h_fluctuation)-1)
GLMakie.scatter!(ax, GLMakie.Point2f0([0.315,1.785]); color=:black, markersize=40, marker='H')
necks2 = Vector{GLMakie.Scatter{Tuple{Vector{Point{2, Float32}}}}}(undef,8)
necks4 = Vector{GLMakie.Scatter{Tuple{Vector{Point{2, Float32}}}}}(undef,8)
necks6 = Vector{GLMakie.Scatter{Tuple{Vector{Point{2, Float32}}}}}(undef,8)
for i in 1:length(isotropy2)
    necks2[i] = GLMakie.scatter!(ax, isotropy2[i], fluctuations2[i]; markersize=30, markerstrokewidth=0, color = cgrad(:greens, 16)[Int(round(16*i/length(fluctuations2)))], markerstrokealpha=0, label = "$(i+1) Necks")
    necks4[i] = GLMakie.scatter!(ax, isotropy4[i], fluctuations4[i]; markersize=30, markerstrokewidth=0, color = cgrad(:reds, 16)[Int(round(16*i/length(fluctuations4)))], markerstrokealpha=0, label = "$(i+1) Necks")
    necks6[i] = GLMakie.scatter!(ax, isotropy6[i], fluctuations6[i]; markersize=30, markerstrokewidth=0, color = cgrad(:blues, 16)[Int(round(16*(i)/length(fluctuations6)))], markerstrokealpha=0, label = "$(i+1) Necks")
end
#GLMakie.Legend(plt[1,2], vcat(necks2, necks4), ["$((i-1)%8+1) Necks in $(div(i-1,8)*2+2) Layers" for i in 1:16])
display(plt)
save("disorderedfluctuations.png", plt)