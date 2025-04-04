using DrWatson, Plots, LaTeXStrings
using Plots.PlotMeasures

gr()  # Use GR backend, but you can switch to PyPlot if needed
function ticks_length!(; tl = 0.02)
	p = Plots.current()
	xticks, yticks = Plots.xticks(p)[1][1], Plots.yticks(p)[1][1]
	xl, yl = Plots.xlims(p), Plots.ylims(p)
	x1, y1 = zero(yticks) .+ xl[1], zero(xticks) .+ yl[1]
	sz = p.attr[:size]
	r = sz[1] / sz[2]
	dx, dy = tl * (xl[2] - xl[1]), tl * r * (yl[2] - yl[1])
	plot!([xticks xticks]', [y1 y1 .+ dy]', c = :black, labels = false)
	plot!([x1 x1 .+ dx]', [yticks yticks]', c = :black, labels = false, xlims = xl, ylims = yl)
	return Plots.current()
end
# Define data


T = 0.1
α = 0.5


# T = 1.0, alpha = 1.5, N_chain = 160, maxdim = 80, cutoff = -11, tau = 0.002, jump = 10, boson_dim = 13, omega = 5
mQ1 = [-1.2250295150933485e-10, 0.0007500588337063434, 0.0896475302875524, 0.31672121848737955, 0.657444843016371, 1.0791427462308694, 1.5478533208036618, 2.0336172890831143, 2.513234938360908, 2.970783805605611, 3.3967334209688835, 3.7864843782528284, 4.138880901143337, 4.454963415882669, 4.737035219680152, 4.9880184112944645, 5.211038300044452, 5.409172692309529, 5.585312890657906, 5.742096661469061, 5.8818854093575945, 6.0067670473095935, 6.118572149775206, 6.218899312500454, 6.309137406428616, 6.3904926080627105, 6.464011359834338, 6.530601650820891, 6.591051941385042, 6.646047791614634, 6.69618632768266, 6.741988753731023, 6.783911100914819, 6.822353435076127, 6.85766714093514, 6.890163767185455, 6.920117518321032, 6.9477720996428225, 6.973344279481128, 6.99702733674857, 7.018993986952362, 7.039398881707598, 7.05838074246547, 7.07606420342883, 7.0925613982716085, 7.107972757431505, 7.122390497483904, 7.135896331865693, 7.148564584379158, 7.160462481035053, 7.171650884206837, 7.182184941889006, 7.19211467124013, 7.201485476558303, 7.210338061426337, 7.218711020594091, 7.226638008416085, 7.234150171019304, 7.2412759304098255, 7.248041253851942, 7.254469895235493, 7.260583610294684, 7.266402352218298, 7.2719438852573095, 7.277226140609878, 7.282264185124312, 7.287072276064765, 7.291663599808316, 7.296050369146868, 7.300244335952117, 7.304255177519679, 7.30809254770433, 7.3117665651984245, 7.315285158935615, 7.318656198560629, 7.321887018617922, 7.324984680293707, 7.327955342204676, 7.330805056621937, 7.333537360405917, 7.336163498706162, 7.338684114522012, 7.341104295091994, 7.343426872056092, 7.345658183110129, 7.34780381828502, 7.3498645659099395, 7.351845117281368, 7.35374795655119, 7.355577941307152, 7.357336889437725, 7.359021440585157, 7.360636387504235, 7.362183875329735, 7.363669254022425, 7.365093497524585, 7.36645913076835, 7.36776883095147, 7.369051678202045, 7.3702643328734645, 7.371427343711353, 7.372541370847603, 7.373613178736233, 7.374687296665454, 7.37569236150371, 7.376669536194154, 7.37761582398932, 7.378523304482236, 7.379389478094884, 7.380220224154634, 7.381045418613671, 7.381795045311904, 7.382499508582544, 7.383151791996622, 7.383768045551927, 7.384356094225934, 7.384922477693131, 7.385494682418238, 7.386039605808739, 7.386601530895642, 7.387618146458499, 7.389389421614149, 7.391803793770689, 7.392378203928341, 7.391074365104866, 7.391061523151933, 7.390648815578756, 7.391708725877975, 7.390210819170091, 7.391711037205523, 7.391604781309504, 7.392441991199415, 7.390928076016276, 7.39122861454506, 7.39159844931756, 7.390799167949074, 7.39101251511068, 7.390973610242075, 7.390121969529086, 7.390344746065562, 7.391938395976837, 7.393420562056871, 7.394072163992797, 7.394687206619547, 7.394987559978177, 7.396052538897253, 7.395071983200764, 7.395656397916173, 7.395787883244343, 7.39459357003488, 7.393399394378623]

vQ1 = [1.4298997150675047e-20, 4.25813647792021e-10, 1.3349391272682527, 4.620739120624743, 9.293364825322291, 14.6484281090442, 20.07000262926544, 25.056616658742826, 29.36615683425323, 32.90952712899926, 35.710687331018384, 37.85451144615432, 39.449054848984986, 40.602687407025115, 41.412546087731045, 41.95861666556884, 42.310363309270656, 42.51711641306956, 42.61855624727882, 42.64475561039209, 42.618282570856515, 42.5548592545562, 42.49158225263809, 42.394843124264206, 42.290323428265836, 42.182724248194624, 42.074631927887275, 41.96996987558962, 41.86918904794095, 41.77324695981969, 41.682049801618625, 41.59732801997297, 41.51817704583109, 41.44452001374464, 41.398497151033816, 41.33745100853495, 41.281737679794915, 41.230567385577174, 41.18363027243938, 41.140402745653, 41.10112456768559, 41.065165887843385, 41.03226304510256, 41.00200696321689, 40.97452965437966, 40.97265850395745, 40.95167386067473, 40.93258137830344, 40.915417320266506, 40.89986049553339, 40.88575180656035, 40.87288470008348, 40.86125767451752, 40.850654732150076, 40.86302818936374, 40.85580921466969, 40.84931809123527, 40.84346368465377, 40.83813456070078, 40.83331394427325, 40.82889041344257, 40.824800615498745, 40.82097523619389, 40.83963078316099, 40.837597675795365, 40.835712962695034, 40.83393044326256, 40.83222592185581, 40.830559817742426, 40.826640986687046, 40.824739276242916, 40.84524983865964, 40.844987418062516, 40.84464409962752, 40.84417266626715, 40.84355759219305, 40.84280586463933, 40.841894389506905, 40.840833495741066, 40.86153270835598, 40.86115411645495, 40.860690576117946, 40.86010001087322, 40.8596296884478, 40.8589364187841, 40.85783213751839, 40.87855309862474, 40.878080922842635, 40.87748219057324, 40.876945128697635, 40.87604781175979, 40.87503623012723, 40.873895124953755, 40.89421077049089, 40.893608875874676, 40.89287732758077, 40.89201154273394, 40.89100615103915, 40.88948159842995, 40.88801545766426, 40.90737466931297, 40.90648447975345, 40.905524969300274, 40.90418280675693, 40.902978822124, 40.90090480128064, 40.898965866063776, 40.896944825088575, 40.91638646494172, 40.91495106648985, 40.91283848911503, 40.9117720627523, 40.90997758280927, 40.90844924532296, 40.90704681173393, 40.92703537996318, 40.92490253895353, 40.92378261104604, 40.92089463352232, 40.91859014120542, 40.98399847509889, 41.29052545198261, 41.55698725809713, 41.59309739336483, 41.24827865979962, 41.065749967399746, 40.902550254343865, 40.993203831919054, 40.670318097217226, 40.78743454150367, 40.709008023609506, 40.870282187026234, 40.71439769426242, 40.75111609371526, 40.88044410444499, 40.729454867108586, 40.71225080619555, 40.71289820847299, 40.60046423870626, 40.59337977015237, 40.83306906188727, 41.00504325011862, 41.02301954134298, 41.141050042375646, 41.20474605111753, 41.45639659528216, 41.37125966178108, 41.49271020888929, 41.36097860731997, 41.083719026354174, 40.83893078567405]


#


#T = 0.1, alpha = 0.1, N_chain = 150, maxdim = 70, cutoff = -11, tau = 0.002, boson_dim = 11;
mQ2 = [1.5155086070002432e-9, 0.00037500713791577164, 0.04482375473748985, 0.1583605863911302, 0.32872245306113657, 0.5395716601368836, 0.7739279326982448, 1.0168114237810397, 1.2566228240897956, 1.4854012991089816, 1.6983820613810479, 1.8932658936639681, 2.0694754073742616, 2.2275312920613164, 2.368585623127957, 2.494099806248687, 2.605636776152643, 2.7047356399653046, 2.7928421842209294, 2.871275357667332, 2.9412158565114632, 3.0037075552320234, 3.0596658575435423, 3.109889312084799, 3.155072312457063, 3.195817650288115, 3.2326482869714743, 3.2660180337129288, 3.2963210648758636, 3.3239002494317282, 3.3490544023147657, 3.3720445253333344, 3.393099160685756, 3.4124189544656676, 3.4301805336305176, 3.4465397972891783, 3.4616347028083734, 3.4755876206342355, 3.488507343106653, 3.5004907817385598, 3.5116244228478086, 3.5219855728383456, 3.531643433348812, 3.5406600278954943, 3.5490910140530247, 3.556986385967239, 3.5643911018007226, 3.571345625644105, 3.5778864153539622, 3.58404635044551, 3.5898551142945436, 3.595339533443954, 3.6005238812265423, 3.6054301445679116, 3.610078266192212, 3.6144863600114783, 3.6186709017285175, 3.622646896543709, 3.6264280336330192, 3.630026824166726, 3.633454716442638, 3.6367222059681756, 3.6398389336923924, 3.64281376697573, 3.6456548784973912, 3.6483698135180385, 3.6509655489803574, 3.653448548902979, 3.6558248114232588, 3.658099913710029, 3.6602790464880077, 3.6623670514825566, 3.664368449860318, 3.6662874716570966, 3.668128079861812, 3.6698939911977884, 3.671588697226915, 3.673215481178506, 3.6747774359829246, 3.676277477621394, 3.6777183566710585, 3.6791026741322788, 3.680432888323256, 3.6817113264265755, 3.6829401945668994, 3.6841215840124684, 3.6852574798534894, 3.686349767892302, 3.6874002393866325, 3.6884105988993374, 3.6893824675228952, 3.6903173918865217, 3.691216841793614, 3.6920822217007836, 3.692914868224074, 3.6937160606766937, 3.6944870167580564, 3.6952289021970404, 3.6959428286468907, 3.6966298610283053, 3.6972910167523216, 3.697927268123742, 3.6985395472941067, 3.6991287452995865, 3.6996957166349262, 3.700241279143616, 3.7007662159063623, 3.70127127798876, 3.7017571857274327, 3.702224630634359, 3.7026742734000795, 3.7031067513735447, 3.7035226733900557, 3.703922625880669, 3.7043071697859373, 3.704676847273795, 3.705032175330946, 3.705373653083167, 3.705701759749896, 3.7060169556935736, 3.7063196834360523, 3.70661036939919, 3.7068894221772335, 3.7071572370319523, 3.7074141923139625, 3.707660653419206, 3.7078969708303204, 3.7081234832073826, 3.7083405167153125, 3.7085483836469666, 3.7087473880978297, 3.7089378199979572, 3.7091199593029094, 3.7092940775390084, 3.7094604344974837, 3.709619280946163, 3.7097708604873003, 3.7099154043533993, 3.710053140006681, 3.7101842848268323, 3.710309047027771, 3.7104276282624924, 3.710540224381721, 3.710647023080712, 3.710748206269503, 3.711593175285782, 3.7135463747758273, 3.7144422325325688, 3.7147616777322066, 3.713618117957227, 3.7122463084724457]
#mQ2 = [-6.56343197184719e-10, 0.00037503152527405045, 0.04482375275418565, 0.15836058447075496, 0.3287224504565828, 0.5395716582732045, 0.7739279302186363, 1.0168114216311068, 1.2566228222414235, 1.4854012970290893, 1.6983820593041112, 1.893265892103666, 2.0694754054813873, 2.2275312897091912, 2.3685856213903427, 2.494099804571196, 2.6056367739956254, 2.704735638359691, 2.7928421820276603, 2.871275355702828, 2.9412158548776897, 3.003707553329315, 3.059665855570393, 3.1098893100891876, 3.1550723103790888, 3.1958176485343444, 3.2326482846703244, 3.2660180324273806, 3.2963210630352515, 3.323900247755791, 3.349054400558976, 3.37204452338141, 3.3930991586218626, 3.4124189521564796, 3.4301805316644387, 3.446539795619061, 3.461634700968578, 3.475587618873872, 3.4885073411372804, 3.500490780133488, 3.5116244211153003, 3.521985571369287, 3.5316434318390804, 3.540660026518644, 3.5490910128506767, 3.556986384798856, 3.5643911003586717, 3.5713456241285852, 3.5778864138265716, 3.5840463489983527, 3.5898551131234018, 3.5953395324924364, 3.600523880406984, 3.605430143421232, 3.6100782652746832, 3.6144863591556966, 3.6186709004914515, 3.6226468947795625, 3.6264280324352223, 3.630026823295787, 3.633454715099589, 3.636722205153093, 3.639838932980315, 3.6428137663112885, 3.645654877874907, 3.648369812258744, 3.6509655479241774, 3.6534485481705214, 3.6558248114592202, 3.658099913478267, 3.660279045999805, 3.6623670510703605, 3.6643684494156696, 3.6662874710961484, 3.6681280795303755, 3.6698939907471044, 3.671588696522402, 3.6732154806413515, 3.674777435804136, 3.6762774771976487, 3.6777183568289304, 3.679102674131592, 3.6804328880023736, 3.6817113263483154, 3.6829401951855734, 3.684121584663594, 3.6852574809630725, 3.6863497684132187, 3.687400240119983, 3.688410599436053, 3.6893824685109804, 3.6903173928065867, 3.691216842708835, 3.6920822222990113, 3.692914869612758, 3.693716062014533, 3.6944870180636107, 3.695228903269036, 3.6959428306233257, 3.696629862734151, 3.6972910183971344, 3.6979272704301573, 3.698539549876959, 3.699128748669191, 3.6996957212816515, 3.700241283835293, 3.7007662221713358, 3.701271285873336, 3.701757196723344, 3.7022246440352076, 3.7026742911138313, 3.7031067727499707, 3.7035227003621447, 3.7039226586614395, 3.7043072107632784, 3.70467689536199, 3.7050322337216484, 3.705373722061584, 3.7057018402302146, 3.7060170483695125, 3.7063197890407276, 3.7066104880580997, 3.7068895547685776, 3.707157382127709, 3.7074143491901066, 3.7076608203250445, 3.707897145680664, 3.7081236632877426, 3.7082266723658317, 3.70953279640737, 3.711476799128999, 3.711379994870493, 3.7105515810448075, 3.708907276563045, 3.707138795402542, 3.7071912270657577, 3.7078046997910232, 3.7078752359923364, 3.7082947171315834, 3.7089338158392, 3.709288187256912, 3.7098015238864286, 3.7104599035683274, 3.709584330025628, 3.708814964045356, 3.708893923821162, 3.7088541671905526, 3.708288165881564, 3.708309473479262, 3.7081629071601605, 3.7078295907740006]
vQ2 = [3.674184411692731e-20, 8.882156347755832e-11, 0.6668659034845279, 2.3082212034471383, 4.641993185942718, 7.316019196107868, 10.00990712315076, 12.492650102942696, 14.635167480693259, 16.39299158177842, 17.778178255666514, 18.833250144920406, 19.61232570244834, 20.169686447208218, 20.552396302748768, 20.805046492484937, 20.958405924344866, 21.03839113328797, 21.064954528730553, 21.051637477500527, 21.01360334738765, 20.956991969933554, 20.888146483681155, 20.811723507052207, 20.730016207052696, 20.64805660928561, 20.565947594265335, 20.484956102433856, 20.404983039804982, 20.328858089005855, 20.255492079387288, 20.18512708422849, 20.11713977914054, 20.05331371895041, 19.992500196165203, 19.934670944682907, 19.87923067871329, 19.827331403891712, 19.778068870416543, 19.73135410641265, 19.686738853297243, 19.6449368096873, 19.605306427154154, 19.56740999695907, 19.531972905926324, 19.498366116417593, 19.466508001535075, 19.436106466402812, 19.40756448615606, 19.380480471752833, 19.35459265895674, 19.33027506867132, 19.30716720057869, 19.285207818070685, 19.26423260511733, 19.244424251396108, 19.22556600330062, 19.207520845002016, 19.190444678714037, 19.17415423152138, 19.15860938156912, 19.1437264889473, 19.129567230546563, 19.11602778765439, 19.10304177254909, 19.09065617962221, 19.078788576499427, 19.06738725219862, 19.05648478943632, 19.04601814066208, 19.035947483279507, 19.026293965303083, 19.017009279890402, 19.00807531830956, 18.999467813647506, 18.991186721017684, 18.983204994296663, 18.97550429378498, 18.96808121579221, 18.960915847513967, 18.953993712493897, 18.94730962962243, 18.94084862103124, 18.93459943444615, 18.928555807447708, 18.922706509667144, 18.91704258044872, 18.91155762558291, 18.90624298848813, 18.90109129443772, 18.896096493290127, 18.89125170142184, 18.886550823952593, 18.881988322685324, 18.877558613716168, 18.873256524576526, 18.869077199392414, 18.86501594014774, 18.86106832297933, 18.857230131914335, 18.853497339907307, 18.84986613382108, 18.846332867506167, 18.842894055063862, 18.83954638655651, 18.83628669362203, 18.833111943727744, 18.830019236029134, 18.827005815138286, 18.824069029442725, 18.82120633843501, 18.8184153095005, 18.81569362292801, 18.81303904296727, 18.810449434648355, 18.807922739391188, 18.80545698203123, 18.8030502750948, 18.800700788232923, 18.798406781553552, 18.796166578723746, 18.79397855363509, 18.79184115891651, 18.78975288939112, 18.78771231133146, 18.785718040418338, 18.78376874348779, 18.781863136120787, 18.77999998157916, 18.778178098358794, 18.776396360387963, 18.774653678498723, 18.772948997475538, 18.771281308065404, 18.769649640223808, 18.768053089717615, 18.766490759549967, 18.764961824131102, 18.76346546644769, 18.762000898083812, 18.76056735503872, 18.759164100799616, 18.75779043463618, 18.756445711696895, 18.755129382356927, 18.87448180389558, 19.1520391644711, 19.253018519441895, 19.20602512298513, 18.946618034423985, 18.71072047628091]
#vQ2 = [9.570858894319533e-20, 1.0840571227784381e-10, 0.6668849505965553, 2.308231445670484, 4.642049265137984, 7.316156210526322, 10.010129692639346, 12.492933694094063, 14.635475729683069, 16.393289632067933, 17.778444958089953, 18.83347475446395, 19.612506427659504, 20.169826942401464, 20.553264126253445, 20.805665211728126, 20.95884122516876, 21.038693463048766, 21.06516380717286, 21.05258384420283, 21.014236000012158, 20.957413177310364, 20.888426002926316, 20.811908653393772, 20.730714385446085, 20.648508663781563, 20.566241566711955, 20.48514696438851, 20.405630234550202, 20.329273589962533, 20.255758562067633, 20.185298055166406, 20.1176716692174, 20.05365080586374, 19.992715096050688, 19.934807804068182, 19.879630755041617, 19.827584860369, 19.778229247333286, 19.73145561175363, 19.687019014655696, 19.645111630786992, 19.605416810896074, 19.567702860561166, 19.532157189103096, 19.498481860812156, 19.466580559550568, 19.436291245850008, 19.407679755535302, 19.380551618926685, 19.35476810928737, 19.3303847747748, 19.307235604815073, 19.28525034000677, 19.264334038124133, 19.244487041152453, 19.22560466564223, 19.20761000694353, 19.19049880501461, 19.174187651925433, 19.15862992792782, 19.1437728689983, 19.129595566984374, 19.11604500370323, 19.103079616984544, 19.090678829378675, 19.078801963490633, 19.067416051818874, 19.056502164131842, 19.046028544469063, 19.03596925911009, 19.026306854270214, 19.017016848250844, 19.00807967491876, 18.999476580959684, 18.991191633499625, 18.983207848884952, 18.9755099649974, 18.968084445885232, 18.960917660483837, 18.953997169622532, 18.947311511478404, 18.9408496225274, 18.934601239704563, 18.92855674183577, 18.922707013920444, 18.9170434879835, 18.911558097083287, 18.906243222940223, 18.90109168405366, 18.896096682783863, 18.89125179991361, 18.886550952479723, 18.88198838779011, 18.87755865459276, 18.873256573066257, 18.86907723440058, 18.865015974486518, 18.861068358341324, 18.857230168869943, 18.853497393209167, 18.849866214355025, 18.846332986728967, 18.842894240409954, 18.839546653830375, 18.836287071300127, 18.833112459440837, 18.83001994281373, 18.827006757798486, 18.82407026661721, 18.821207941792167, 18.81841735967234, 18.815696207865958, 18.81304227862717, 18.81045345273365, 18.8079277072705, 18.805463108608677, 18.80305782566398, 18.800710089120738, 18.798418229484053, 18.796180654611042, 18.793995852323775, 18.791862357167687, 18.789778780589888, 18.78774380731702, 18.785756141382443, 18.78381457663319, 18.78191798139563, 18.76168314869348, 18.94352620088346, 19.216212056932505, 19.137918194565575, 19.001473611084705, 18.783563827818202, 18.620198800769142, 18.585193752405537, 18.71554452778911, 18.679818143551437, 18.547665394552325, 18.58483622906374, 18.711770820564027, 18.87536509174428, 18.941123622700182, 18.774392296937826, 18.658590997710515, 18.730667368557192, 18.733547275777546, 18.66340768219719, 18.692257238480465, 18.66386293245346, 18.679971441003126]

#T = 1, alpha = 0.75, N_chain = 160, maxdim = 80, cutoff = -11, tau = 0.002, jump = 10, boson_dim = 14, omega = 5
mQ3 = [-1.5734394870880753e-10, 0.00037502139792946316, 0.04482376429483656, 0.15836060616513828, 0.3287224150148531, 0.5395713620062367, 0.7739269563738549, 1.0168090727018446, 1.2566180764840855, 1.4853920607818982, 1.698367160186902, 1.8932432588280996, 2.0694426595536233, 2.227485782410144, 2.3685244744524363, 2.494019954338646, 2.605535007635483, 2.7046086290408367, 2.792686519019119, 2.871087573641152, 2.94099246592658, 3.003445080657642, 3.05936086448405, 3.109538435730957, 3.1546722725901035, 3.195365249740921, 3.2321403864538683, 3.2654515093704264, 3.2956927386521366, 3.323206212901304, 3.3482916837673846, 3.371209318782723, 3.3921873071618567, 3.411425898038437, 3.4291012985157012, 3.4453689932531573, 3.4603665517773847, 3.4742160078557043, 3.487025884104446, 3.4988929052295275, 3.5099034610388062, 3.5201348552992746, 3.5296563810367365, 3.5385302438774042, 3.546811748843059, 3.5545524240775355, 3.5617970531347853, 3.5685865555212506, 3.574957887692708, 3.5809444565772064, 3.586576493798658, 3.5918813849619604, 3.5968839614579156, 3.6016067651208363, 3.606070279189978, 3.6102931402755485, 3.6142923228796113, 3.618082711955473, 3.6216796035130843, 3.625095352732267, 3.6283417894317496, 3.631429760034772, 3.634369218861515, 3.6371693177545126, 3.639838481080982, 3.6423844747659357, 3.6448144682235033, 3.6471350888350416, 3.6493524726095568, 3.651472306074996, 3.6534992862315216, 3.655439471064102, 3.657296862531289, 3.659075718668298, 3.66078001333787, 3.662413458444197, 3.6639795255635046, 3.6654814674731613, 3.6669223332315326, 3.668304984561657, 3.669632109855978, 3.670906239320762, 3.672129752725333, 3.6733048947163245, 3.674433195208073, 3.6755178122401215, 3.6765600599011736, 3.6775617225796013, 3.678524489445574, 3.679449960682206, 3.6803396507107613, 3.6811949977415326, 3.6820173663249425, 3.6828080510387466, 3.68356828340833, 3.684299232914364, 3.6850020131214976, 3.6856776837212872, 3.68632667201778, 3.68695109341923, 3.687551291584923, 3.6881281404924184, 3.6886824734679884, 3.68921508696132, 3.6897267411323176, 3.6902181619716186, 3.6906900441002555, 3.6911430510000494, 3.691577818010488, 3.6919949532194707, 3.6923950378444834, 3.692778629081224, 3.6931456884844014, 3.6934978691004243, 3.693835092747026, 3.694157830339096, 3.6944665326237143, 3.6947616335781457, 3.695043549785443, 3.695312677733939, 3.6955694036297824, 3.695814095489376, 3.6960471073159975, 3.6962687798533462, 3.696479440649476, 3.696679403735934, 3.696868973412044, 3.6970484402808537, 3.6973664148113663, 3.6961819130847373, 3.694761690075027, 3.6951843409930554, 3.695373437707646, 3.6978211727512327, 3.701318202707647, 3.701276697606342, 3.700930746372221, 3.7011305670729575, 3.700296766213168, 3.69910864533066, 3.700222738394754, 3.7018456421965618, 3.7021590057396456, 3.7018544795626287, 3.702751384987736, 3.703512517546288, 3.7039016212185563, 3.704195987007214, 3.703880265726499, 3.7036878022160344, 3.7040135207235982]
vQ3 = [9.255987564464473e-21, 1.0056776464157568e-10, 0.6674666522547563, 2.3103695321180835, 4.646685090055164, 7.324223547098463, 10.022551783901571, 12.510616610563284, 14.6592973888012, 16.454801931928454, 17.855397823689028, 18.92733670354271, 19.72465068268411, 20.301528067441296, 20.705527848441527, 20.97978289362942, 21.155760426368307, 21.259262571637148, 21.31012520696281, 21.323385443657337, 21.309696893237316, 21.278908950298728, 21.23627770464821, 21.186347635498652, 21.13179761815546, 21.076431161713895, 21.02097645198699, 20.96659855436532, 20.913563033960536, 20.889064058229586, 20.843997271091087, 20.80168633108216, 20.76179538265984, 20.725165325391547, 20.69113430618836, 20.659606457572295, 20.630206941966225, 20.603375438283948, 20.5785916884656, 20.55572242513187, 20.53448225313754, 20.515098333690368, 20.49721006885169, 20.480704585188295, 20.491410655619035, 20.478992921989434, 20.467627383267693, 20.457143756473506, 20.447636252179166, 20.43890181687952, 20.430867150456045, 20.423424964212302, 20.4166014251509, 20.410277330187558, 20.404369141536478, 20.39890190127307, 20.393782037657044, 20.41402932706928, 20.410616037186244, 20.407419782422906, 20.40441016253805, 20.401549084570178, 20.398834130859143, 20.396226772681167, 20.393700887891463, 20.391254496031436, 20.388861519677295, 20.386505033206838, 20.38418318207064, 20.381879631409802, 20.404249817125546, 20.40272352903368, 20.40117500455406, 20.399596274337657, 20.39798204325252, 20.39632991903203, 20.394636055862577, 20.39289813914813, 20.391115366034605, 20.389286514081284, 20.38741130132947, 20.385490146375716, 20.38352353198424, 20.381512463349143, 20.404159068769946, 20.402557590207813, 20.40089648156666, 20.399177442146367, 20.397402349350486, 20.395573278690943, 20.393692433649335, 20.39176213851845, 20.389784811856718, 20.387762927836683, 20.385699007713047, 20.38359559867675, 20.38145525312833, 20.379280511385446, 20.401521234156455, 20.399594524159166, 20.397625122579925, 20.395615582813935, 20.393568437698395, 20.391486169015227, 20.38937124779062, 20.387226083207015, 20.385053032864192, 20.382854407834337, 20.38063246413526, 20.378389388527463, 20.37612730862383, 20.37384828226284, 20.39565339691172, 20.393534017554103, 20.391391546191887, 20.38922801134691, 20.38704537247253, 20.38484554828881, 20.382630397627825, 20.380401717101474, 20.3781612443617, 20.375910702642912, 20.3736517384053, 20.37138597297138, 20.369115030486025, 20.366840496128162, 20.3645639447924, 20.3622869427026, 20.383924777017363, 20.152480866517923, 19.91674309830129, 20.013842328419457, 20.00123147083866, 20.28620909896768, 20.7121556907827, 20.619630616005427, 20.48968187459232, 20.479466016135323, 20.42003778185485, 20.1798444469608, 20.19611276769669, 20.349048399120377, 20.40672409945664, 20.272247252312482, 20.434323469849282, 20.57687137291143, 20.65873860382073, 20.682101489450226, 20.583169859171758, 20.502830613515115, 20.560667604656857]

ω_C = 5
f(t, α) = α*ω_C*(t)^2 /(1+t^2)
#= 
g(ω, t, β) = (1-cos(ω*t)*coth(β*ω/2)*2*α*ω*exp(-ω/ω_C))
vQ(t, α, β) = 0.5 .* [quadgk(g(ω, t[i], β), 0, 10^7) for i in 1:lastindex(t)] =#
#= 
g(ω, t, β, α) = (1 - cos(ω * t)) * coth(β * ω / 2) * 2 * α * ω * exp(-ω / ω_C)
vQ(t, α, β) = 0.5 .* [quadgk(ω -> g(ω, t[i], β, α), 0, 10^7)[2] for i in 1:lastindex(t)]
 =#


function vQ(t, α, β, ω_C=5)
    g(ω, t) = (1 - cos(ω * t/ω_C)) * coth(β * ω / 2) * 2 * α * ω * exp(-ω / ω_C)
    return 0.5 .* [quadgk(ω -> g(ω, t[i]), 0, 1e5)[1] for i in eachindex(t)]
end


# Time step duration
tau = 0.002  
nt = 500
ttotal = nt * tau  # Total time evolution
#time_steps = collect(0:10*tau:ttotal)

time_steps = collect(0:lastindex(mQ3)-1)*10*tau*ω_C

mQ1_exact = f.(time_steps, 1.5)
mQ2_exact = f.(time_steps, 0.75)
mQ3_exact = f.(time_steps, 0.75)

vQ1_exact = vQ(time_steps, 1.5, 1)
vQ2_exact = vQ(time_steps, .75, 10)
vQ3_exact = vQ(time_steps, .75, 1)


# Set tick positions and convert tick labels to LaTeX strings automatically
xticks1 = range(2, stop = 8, length = 2)
xticks2 = range(2, stop = 8, length = 2)
xticks3 = range(2, stop = maximum(time_steps[end-30]), length = 3)
xticks3sub = range(8, stop = maximum(time_steps[end-40]), length = 2)

yticks1 = range(0, stop = 7, length = 3)
yticks2 = range(0, stop = 40, length = 3)
yticks3 = range(5, stop = 8, length = 2)
yticks3sub = range(5, stop = 5.7, length = 2)

#xtick_labels1 = [string(round(x, digits = 1)) for x in xticks1]
xtick_labels1 = [" " for x in xticks1]

xtick_labels2 = [string(round(x, digits = 1)) for x in xticks2]
xtick_labels3 = [string((round(x, digits = 1))) for x in xticks3]
xtick_labels3sub = [string((round(x, digits = 1))) for x in xticks3sub]
ytick_labels1 = [string((round(y, digits = 1))) for y in yticks1]
ytick_labels2 = [string((round(y, digits = 1))) for y in yticks2]
ytick_labels3 = [string((round(y, digits = 1))) for y in yticks3]
ytick_labels3sub = [string((round(y, digits = 1))) for y in yticks3sub]


# Plot
p1 = plot(time_steps, mQ1_exact, color = :teal, seriesalpha = 1, linewidth = 2, label = false)
plot!(time_steps, mQ3_exact, color = :lightcoral, seriesalpha = 1, linewidth = 2, label = false)
plot!(time_steps, mQ2_exact, color = :lightblue, seriesalpha = 1, linewidth = 2, label = false)
scatter!(time_steps[1:4:end], mQ1[1:4:end], markersize = 4, markerstrokewidth = 0.5, xaxis = "", #yaxis = L"\langle Q \rangle", 
	xticks = (xticks1, xtick_labels1), yticks = (yticks1, ytick_labels1),
	xtickfont = font(14), ytickfont = font(14),
	xguidefontsize = 24, yguidefontsize = 22, color = :teal, label = L"\  (1, 1.5)") 
scatter!(time_steps[1:6:end], mQ3[1:6:end], markersize = 4, markerstrokewidth = 0.5, color = :lightcoral,label = L"\ (1, 0.75)")
scatter!(time_steps[4:6:end], mQ2[4:6:end], markersize = 4, markerstrokewidth = 0.5, color = :lightblue, label = L"\ (0.1, 0.75)")
plot!(yaxis = L"⟨Q⟩")
plot!(legend = false)
xlims!(0, 10)  
ylims!(0, 8)  
plot!(widen = false)
vline!([xlims(p1)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p1)[2]], lc = :black, lw = 2, label = false)
plot!(grid = false)
plot!(bottom_margin = -8mm)
annotate!(9.5, 6, text("(a)", 17, :black, :right))
ticks_length!(tl=.03)




p2 = plot(time_steps, vQ1_exact, color = :teal, seriesalpha = 1, linewidth = 2, label = false)
plot!(time_steps, vQ3_exact, color = :lightcoral, seriesalpha = 1, linewidth = 2, label = false)
plot!(time_steps, vQ2_exact, color = :lightblue, seriesalpha = 1, linewidth = 2, label = false)
scatter!(time_steps[1:3:end], vQ1[1:3:end], markersize = 4, markerstrokewidth = 0.5, xticks = (xticks2, xtick_labels2), yticks = (yticks2, ytick_labels2), xtickfont = font(14), ytickfont = font(14),
	xguidefontsize = 24, yguidefontsize = 22, color = :teal, seriesalpha = 1, xaxis = "") #, label = L"\  (1, 1.5)")
scatter!(time_steps[4:5:end], vQ3[4:5:end], markersize = 4, markerstrokewidth = 0.5, color = :lightcoral) #, label = L"\ (0.1, 0.75)")
scatter!(time_steps[2:5:end], vQ2[2:5:end], markersize = 4, markerstrokewidth = 0.5, color = :lightblue)
plot!(widen = false)
xlims!(0, 10)  
ylims!(0, 45)  
plot!(yaxis = L"⟨⟨Q^2⟩⟩")
plot!(xaxislabel = false)
vline!([xlims(p2)[2]], lc = :black, lw = 2, label = false)
hline!([ylims(p2)[2]], lc = :black, lw = 2, label = false)
plot!(grid = false)
plot!(legend = false)
annotate!(9.5, 35, text("(b)", 17, :black, :right))
ticks_length!(tl=.03)
plot!(xaxis = L"tω_C")



p3 = plot(time_steps[10:end], vQ1_exact[10:end]./mQ1_exact[10:end], color = :teal, seriesalpha = 1, linewidth = 2, label = false)
plot!(time_steps[10:end], vQ3_exact[10:end]./mQ3_exact[10:end], color = :lightcoral, seriesalpha = 1, linewidth = 2, label = false)
plot!(time_steps[10:end], vQ2_exact[10:end]./mQ2_exact[10:end], color = :lightblue, seriesalpha = 1, linewidth = 2, label = false)
scatter!(time_steps[8:6:length(vQ1)], vQ1[8:6:end] ./ mQ1[8:6:end], color = :teal, xticks = (xticks3, xtick_labels3), yticks = (yticks3, ytick_labels3), xtickfont = font(14), ytickfont = font(14),
	xguidefontsize = 24, yguidefontsize = 22, seriesalpha = 1, xaxis = "", markersize = 4, markerstrokewidth = 0.5) #, label = L"\  (1, 1.5)")
scatter!(time_steps[8:6:length(vQ2)], vQ2[8:6:end] ./ mQ2[8:6:end], markersize = 4, markerstrokewidth = 0.5, color = :lightblue) #, label = L"\ (0.1, 0.75)")
scatter!(time_steps[8:6:length(vQ3)], vQ3[8:6:end] ./ mQ3[8:6:end], markersize = 4, markerstrokewidth = 0.5, color = :lightcoral)
plot!(widen = false)
#plot!(yaxis = L"F")
ylims!(4.5, 10)  
xlims!(0, 13)  
vline!([xlims(p3)[1]], lc = :black, lw = 2, label = false)
hline!([ylims(p3)[2]], lc = :black, lw = 2, label = false)
plot!(grid = false, ymirror = true)
plot!(legend = false)
annotate!(15, 7, text(L"F", 24, :black, :right))
ticks_length!(tl=.02)
plot!(xaxis = L"tω_C")
annotate!(11, 9.7 , text("(c)", 17, :black, :right)) 


plot!(p3, inset=bbox(0.33,0.06,0.4, 0.4), subplot=2)

plot!(p3[2], time_steps[50:end], vQ1_exact[50:end]./mQ1_exact[50:end], color = :teal, seriesalpha = 1, linewidth = 2, label = false)
plot!(p3[2], time_steps[50:end], vQ3_exact[50:end]./mQ3_exact[50:end], color = :lightcoral, seriesalpha = 1, linewidth = 2, label = false)
plot!(p3[2], time_steps[50:end], vQ2_exact[50:end]./mQ2_exact[50:end], color = :lightblue, seriesalpha = 1, linewidth = 2, label = false)
scatter!(p3[2], time_steps[50:6:length(vQ1)], vQ1[50:6:end] ./ mQ1[50:6:end], color = :teal, xticks = (xticks3sub, xtick_labels3sub), yticks = (yticks3sub, ytick_labels3sub), xtickfont = font(12), ytickfont = font(12),
	xguidefontsize = 10, yguidefontsize = 15, seriesalpha = 1, xaxis = "", markersize = 4, markerstrokewidth = 0.5) #, label = L"\  (1, 1.5)")
scatter!(p3[2], time_steps[50:6:length(vQ3)], vQ3[50:6:end] ./ mQ3[50:6:end], markersize = 4, markerstrokewidth = 0.5, color = :lightcoral) #, label = L"\ (0.1, 0.75)")
scatter!(p3[2], time_steps[50:6:length(vQ2)], vQ2[50:6:end] ./ mQ2[50:6:end], markersize = 4, markerstrokewidth = 0.5, color = :lightblue)
ylims!(p3[2], 5, 5.8)  
xlims!(p3[2], 7, 12)  
#= vline!([xlims(p3)[1]], lc = :black, lw = 2, label = false)
hline!([ylims(p3)[2]], lc = :black, lw = 2, label = false) =#
plot!(p3[2], grid = false)
plot!(p3[2], legend = false, aspect_ratio = 5., framestyle = :box)	
#annotate!(19.7, 4.7, text(L"F", 24, :black, :right))
#annotate!(14.5, 5.8, text("(c)", 20, :black, :right))
#plot!(p3[2], xaxis = L"tω_C")
#plot!(p3[2], xaxis = " ")



custom_layout = @layout [[a{0.5h}; b{1.35w}] c{0.6w}]
p = plot(p1, p2, p3, layout = custom_layout, size = (600, 400), left_margin = 5mm, right_margin = 7mm)



#p=plot(p1, p2, p3, p4, layout=(2,2))

display("image/png", p)

