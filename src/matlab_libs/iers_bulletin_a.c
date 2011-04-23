#include <stdio.h>
#include <assert.h>
#include "iers_bulletin_a.h"

typedef struct { float x; float y; float dut1; } xyz;

const static xyz s1[] = {
  {.197846, .261896, -.1656442}, // 54344
  {.196690, .259844, -.1664020}, // 54345
  {.195332, .257834, -.1668913}, // 54346
  {.194017, .255871, -.1671245}, // 54347
  {.192566, .253884, -.1671878}, // 54348
  {.190582, .251822, -.1672302}, // 54349
  {.188552, .249788, -.1673569}, // 54350
  {.186850, .247717, -.1676501}, // 54351
  {.185549, .246031, -.1681548}, // 54352
  {.184318, .244405, -.1688723}, // 54353
  {.182768, .242604, -.1697019}, // 54354
  {.181114, .240563, -.1705973}, // 54355
  {.179227, .238421, -.1715590}, // 54356
  {.177308, .236307, -.1724556}, // 54357
  {.175077, .234251, -.1732340}, // 54358
  {.172349, .232232, -.1738688}, // 54359
  {.169627, .230326, -.1743528}, // 54360
  {.167222, .228494, -.1747153}, // 54361
  {.165159, .226675, -.1749653}, // 54362
  {.163088, .224793, -.1751827}, // 54363
  {.160462, .223015, -.1754827}, // 54364
  {.157620, .221118, -.1759385}, // 54365
  {.154711, .219589, -.1766034}, // 54366
  {.151735, .217998, -.1775500}, // 54367
  {.149034, .216380, -.1787957}, // 54368
  {.146531, .214823, -.1803971}, // 54369
  {.144031, .213052, -.1821224}, // 54370
  {.141497, .211101, -.1837554}, // 54371
  {.138683, .209497, -.1851753}, // 54372
  {.136108, .207821, -.1863046}, // 54373
  {.133566, .206447, -.1871348}, // 54374
  {.130586, .205313, -.1877198}, // 54375
  {.127516, .204009, -.1882266}, // 54376
  {.124523, .202867, -.1887602}, // 54377
  {.121546, .201685, -.1893749}, // 54378
  {.118595, .200523, -.1900910}, // 54379
  {.115681, .199345, -.1909531}, // 54380
  {.112797, .198279, -.1919368}, // 54381
  {.109778, .197525, -.1929827}, // 54382
  {.106615, .196831, -.1940630}, // 54383
  {.103296, .195956, -.1951621}, // 54384
  {.100232, .195044, -.1961659}, // 54385
  {.097521, .194272, -.1970283}, // 54386
  {.094852, .193584, -.1977069}, // 54387
  {.092006, .192992, -.1982271}, // 54388
  {.088927, .192567, -.1986529}, // 54389
  {.085841, .192103, -.1990366}, // 54390
  {.082942, .191867, -.1994848}, // 54391
  {.079904, .191405, -.2000408}, // 54392
  {.076608, .191030, -.2007568}, // 54393
  {.073424, .190526, -.2016848}, // 54394
  {.070750, .190159, -.2028326}, // 54395
  {.068560, .190154, -.2042226}, // 54396
  {.066393, .190392, -.2058425}, // 54397
  {.064116, .190404, -.2075348}, // 54398
  {.061512, .190078, -.2090663}, // 54399
  {.058856, .190160, -.2102561}, // 54400
  {.056182, .190273, -.2110852}, // 54401
  {.053428, .190423, -.2116676}, // 54402
  {.050650, .190668, -.2121596}, // 54403
  {.047663, .190909, -.2126471}, // 54404
  {.044814, .191049, -.2132040}, // 54405
  {.042285, .191108, -.2138903}, // 54406
  {.040288, .191493, -.2147402}, // 54407
  {.038221, .192127, -.2157262}, // 54408
  {.035920, .192552, -.2167481}, // 54409
  {.034039, .192813, -.2177933}, // 54410
  {.032254, .193090, -.2188466}, // 54411
  {.030184, .193263, -.2199019}, // 54412
  {.027867, .193383, -.2208456}, // 54413
  {.025274, .193441, -.2216030}, // 54414
  {.022748, .193570, -.2221761}, // 54415
  {.019964, .194153, -.2225899}, // 54416
  {.016908, .194712, -.2229131}, // 54417
  {.014037, .194955, -.2232638}, // 54418
  {.011073, .195484, -.2237320}, // 54419
  {.007872, .196064, -.2244439}, // 54420
  {.004561, .196430, -.2254330}, // 54421
  {.001372, .197029, -.2267030}, // 54422
  {-.001555, .198124, -.2282814}, // 54423
  {-.004213, .199402, -.2301190}, // 54424
  {-.006928, .200802, -.2320994}, // 54425
  {-.009466, .202306, -.2340995}, // 54426
  {-.011938, .203726, -.2359422}, // 54427
  {-.014899, .204513, -.2375090}, // 54428
  {-.018139, .204854, -.2387677}, // 54429
  {-.020843, .205170, -.2397811}, // 54430
  {-.023072, .205894, -.2406888}, // 54431
  {-.025300, .206944, -.2416482}, // 54432
  {-.027277, .208180, -.2427560}, // 54433
  {-.028727, .209445, -.2440641}, // 54434
  {-.029633, .211083, -.2455588}, // 54435
  {-.030691, .212637, -.2471393}, // 54436
  {-.032054, .214201, -.2487433}, // 54437
  {-.034042, .215879, -.2502569}, // 54438
  {-.036753, .217196, -.2515918}, // 54439
  {-.039716, .218573, -.2526780}, // 54440
  {-.042368, .219829, -.2535109}, // 54441
  {-.044568, .221124, -.2540921}, // 54442
  {-.046272, .222438, -.2544507}, // 54443
  {-.047708, .223984, -.2546589}, // 54444
  {-.049165, .225523, -.2548138}, // 54445
  {-.050320, .226681, -.2550384}, // 54446
  {-.051632, .227770, -.2553954}, // 54447
  {-.053206, .228809, -.2559116}, // 54448
  {-.054691, .230280, -.2566430}, // 54449
  {-.056181, .231826, -.2576299}, // 54450
  {-.057584, .233411, -.2588381}, // 54451
  {-.058957, .235076, -.2601591}, // 54452
  {-.060065, .236346, -.2614713}, // 54453
  {-.061021, .237712, -.2627366}, // 54454
  {-.061936, .239355, -.2638367}, // 54455
  {-.062800, .241016, -.2646203}, // 54456
  {-.063869, .242720, -.2651745}, // 54457
  {-.065327, .244486, -.2656467}, // 54458
  {-.067349, .246143, -.2661617}, // 54459
  {-.069530, .247660, -.2668088}, // 54460
  {-.071332, .249063, -.2676291}, // 54461
  {-.073061, .251084, -.2686237}, // 54462
  {-.074696, .252765, -.2697643}, // 54463
  {-.076492, .254818, -.2710246}, // 54464
  {-.078584, .256878, -.2722461}, // 54465
};

const int s1_mjd_range[2] = { 54344, 54465 };

const static xyz s2[] = {
  {.119926, .539795, -.4287334}, // 54618
  {.122599, .539282, -.4296733}, // 54619
  {.125604, .538545, -.4303358}, // 54620
  {.129022, .537432, -.4307708}, // 54621
  {.132770, .536157, -.4310959}, // 54622
  {.136409, .535023, -.4314450}, // 54623
  {.139862, .533821, -.4319064}, // 54624
  {.143211, .532810, -.4325315}, // 54625
  {.146566, .531717, -.4333207}, // 54626
  {.150059, .530750, -.4342625}, // 54627
  {.153708, .529230, -.4352712}, // 54628
  {.157068, .527722, -.4361896}, // 54629
  {.160557, .526139, -.4369375}, // 54630
  {.163853, .524745, -.4374952}, // 54631
  {.166759, .523375, -.4378270}, // 54632
  {.169548, .521967, -.4379951}, // 54633
  {.172402, .520748, -.4380640}, // 54634
  {.175154, .519104, -.4380358}, // 54635
  {.177759, .517479, -.4379648}, // 54636
  {.180363, .515939, -.4379085}, // 54637
  {.183588, .514407, -.4379158}, // 54638
  {.186973, .512933, -.4381027}, // 54639
  {.190223, .511392, -.4385269}, // 54640
  {.193525, .510010, -.4392002}, // 54641
  {.196463, .508450, -.4400651}, // 54642
  {.199396, .507144, -.4410398}, // 54643
  {.201886, .505901, -.4420680}, // 54644
  {.204186, .504066, -.4429893}, // 54645
  {.206511, .502126, -.4437771}, // 54646
  {.208734, .499884, -.4443840}, // 54647
  {.211255, .497632, -.4447553}, // 54648
  {.213559, .495370, -.4449413}, // 54649
  {.216149, .493286, -.4450863}, // 54650
  {.218818, .491229, -.4453621}, // 54651
  {.221700, .489137, -.4458309}, // 54652
  {.224624, .487291, -.4464747}, // 54653
  {.227207, .485447, -.4472565}, // 54654
  {.229771, .483384, -.4480664}, // 54655
  {.231853, .481173, -.4487570}, // 54656
  {.234144, .479243, -.4493163}, // 54657
  {.236459, .477280, -.4497288}, // 54658
  {.238787, .475010, -.4498757}, // 54659
  {.241179, .472528, -.4497844}, // 54660
  {.243325, .469909, -.4495398}, // 54661
  {.245291, .467105, -.4491896}, // 54662
  {.247842, .464123, -.4487776}, // 54663
  {.250546, .461605, -.4484082}, // 54664
  {.252916, .459118, -.4481648}, // 54665
  {.254954, .456614, -.4480718}, // 54666
  {.257172, .453836, -.4482027}, // 54667
  {.259583, .450902, -.4486146}, // 54668
  {.262270, .447726, -.4492319}, // 54669
  {.265041, .444176, -.4499528}, // 54670
  {.267752, .440735, -.4507245}, // 54671
  {.270205, .437449, -.4514756}, // 54672
  {.272792, .434829, -.4521019}, // 54673
  {.274957, .432400, -.4525169}, // 54674
  {.276706, .429685, -.4527193}, // 54675
  {.278295, .426967, -.4527730}, // 54676
  {.279393, .424095, -.4527443}, // 54677
  {.280189, .421356, -.4527960}, // 54678
  {.281015, .418309, -.4530106}, // 54679
  {.282233, .415093, -.4534397}, // 54680
  {.283650, .411639, -.4540604}, // 54681
  {.285242, .408117, -.4548430}, // 54682
  {.286912, .404654, -.4556345}, // 54683
  {.288545, .401369, -.4563016}, // 54684
  {.289915, .398428, -.4568091}, // 54685
  {.290768, .395594, -.4571338}, // 54686
  {.291378, .392687, -.4572796}, // 54687
  {.291927, .389503, -.4572741}, // 54688
  {.292402, .386011, -.4571430}, // 54689
  {.292772, .382648, -.4569244}, // 54690
  {.292989, .379366, -.4567038}, // 54691
  {.293028, .376324, -.4565913}, // 54692
  {.293033, .373400, -.4566224}, // 54693
  {.293585, .370602, -.4568064}, // 54694
  {.294201, .367988, -.4571466}, // 54695
  {.294489, .365067, -.4577060}, // 54696
  {.294813, .361841, -.4584595}, // 54697
  {.295375, .358674, -.4592921}, // 54698
  {.296075, .355608, -.4600717}, // 54699
  {0.296453, 0.352475, -.4606745}, // 54700
  {0.296801, 0.349322, -.4610545}, // 54701
  {0.296945, 0.346068, -.4612418}, // 54702
  {0.296939, 0.342717, -.4613117}, // 54703
  {0.296831, 0.339286, -.4613523}, // 54704
  {0.296649, 0.335798, -.4614655}, // 54705
  {0.296402, 0.332275, -.4617495}, // 54706
  {0.296090, 0.328735, -.4622745}, // 54707
  {0.295714, 0.325186, -0.4630471}, // 54708
  {0.295272, 0.321636, -0.4640229}, // 54709
  {0.294765, 0.318086, -0.4650688}, // 54710
  {0.294194, 0.314538, -0.4660528}, // 54711
  {0.293560, 0.310994, -0.4669029}, // 54712
  {0.292864, 0.307456, -0.4675628}, // 54713
  {0.292107, 0.303926, -0.4680198}, // 54714
  {0.291290, 0.300404, -0.4683243}, // 54715
  {0.290413, 0.296892, -0.4685212}, // 54716
  {0.289478, 0.293393, -0.4686735}, // 54717
  {0.288484, 0.289906, -0.4688471}, // 54718
  {0.287431, 0.286432, -0.4691097}, // 54719
  {0.286320, 0.282974, -0.4695332}, // 54720
  {0.285151, 0.279531, -0.4701730}, // 54721
  {0.283925, 0.276105, -0.4710708}, // 54722
  {0.282641, 0.272696, -0.4722353}, // 54723
  {0.281301, 0.269306, -0.4736246}, // 54724
  {0.279905, 0.265935, -0.4751437}, // 54725
  {0.278452, 0.262584, -0.4766613}, // 54726
  {0.276944, 0.259254, -0.4780526}, // 54727
  {0.275380, 0.255946, -0.4792253}, // 54728
  {0.273761, 0.252661, -0.4801554}, // 54729
  {0.272088, 0.249400, -0.4808904}, // 54730
  {0.270361, 0.246162, -0.4815341}, // 54731
  {0.268579, 0.242951, -0.4822169}, // 54732
  {0.266745, 0.239765, -0.4830548}, // 54733
  {0.264858, 0.236606, -0.4841229}, // 54734
  {0.262918, 0.233475, -0.4854393}, // 54735
  {0.260927, 0.230372, -0.4869649}, // 54736
  {0.258884, 0.227298, -0.4886160}, // 54737
  {0.256790, 0.224255, -0.4902879}, // 54738
  {0.254646, 0.221243, -0.4918743}, // 54739
  {0.252452, 0.218262, -0.4933005}, // 54740
  {0.250208, 0.215313, -0.4945303}, // 54741
  {0.247916, 0.212398, -0.4955634}, // 54742
  {0.245576, 0.209517, -0.4964362}, // 54743
  {0.243188, 0.206670, -0.4972076}, // 54744
  {0.240754, 0.203858, -0.4979441}, // 54745
  {0.238272, 0.201083, -0.4987193}, // 54746
  {0.235745, 0.198345, -0.4996058}, // 54747
  {0.233173, 0.195644, -0.5006623}, // 54748
  {0.230557, 0.192981, -0.5019310}, // 54749
  {0.227897, 0.190357, -0.5034337}, // 54750
  {0.225193, 0.187773, -0.5051530}, // 54751
  {0.222447, 0.185229, -0.5070241}, // 54752
  {0.219659, 0.182725, -0.5089366}, // 54753
  {0.216830, 0.180264, -0.5107601}, // 54754
  {0.213961, 0.177844, -0.5123825}, // 54755
  {0.211052, 0.175467, -0.5137450}, // 54756
  {0.208104, 0.173134, -0.5148686}, // 54757
  {0.205118, 0.170845, -0.5158506}, // 54758
  {0.202094, 0.168600, -0.5168249}, // 54759
  {0.199034, 0.166400, -0.5179178}, // 54760
  {0.195938, 0.164246, -0.5192123}, // 54761
  {0.192807, 0.162139, -0.5207290}, // 54762
  {0.189641, 0.160078, -0.5224289}, // 54763
  {0.186442, 0.158065, -0.5242362}, // 54764
  {0.183210, 0.156099, -0.5260581}, // 54765
  {0.179946, 0.154182, -0.5278037}, // 54766
  {0.176650, 0.152314, -0.5293959}, // 54767
  {0.173325, 0.150495, -0.5307899}, // 54768
  {0.169970, 0.148726, -0.5319748}, // 54769
  {0.166587, 0.147008, -0.5329723}, // 54770
  {0.163176, 0.145340, -0.5338323}, // 54771
  {0.159738, 0.143723, -0.5346197}, // 54772
  {0.156274, 0.142158, -0.5354052}, // 54773
  {0.152784, 0.140645, -0.5362565}, // 54774
  {0.149271, 0.139185, -0.5372321}, // 54775
  {0.145734, 0.137777, -0.5383752}, // 54776
  {0.142175, 0.136422, -0.5397090}, // 54777
  {0.138594, 0.135121, -0.5412344}, // 54778
  {0.134992, 0.133874, -0.5429174}, // 54779
  {0.131371, 0.132681, -0.5446834}, // 54780
  {0.127731, 0.131543, -0.5464202}, // 54781
  {0.124073, 0.130459, -0.5480074}, // 54782
  {0.120397, 0.129430, -0.5493565}, // 54783
  {0.116706, 0.128457, -0.5504457}, // 54784
  {0.113000, 0.127539, -0.5513367}, // 54785
  {0.109279, 0.126678, -0.5521535}, // 54786
  {0.105545, 0.125872, -0.5530342}, // 54787
  {0.101799, 0.125123, -0.5540809}, // 54788
  {0.098041, 0.124430, -0.5553304}, // 54789
  {0.094273, 0.123794, -0.5567506}, // 54790
  {0.090496, 0.123214, -0.5582629}, // 54791
  {0.086710, 0.122692, -0.5597717}, // 54792
  {0.082917, 0.122227, -0.5611889}, // 54793
  {0.079117, 0.121819, -0.5624478}, // 54794
  {0.075311, 0.121469, -0.5635084}, // 54795
  {0.071500, 0.121176, -0.5643584}, // 54796
  {0.067686, 0.120941, -0.5650124}, // 54797
  {0.063870, 0.120763, -0.5655111}, // 54798
  {0.060051, 0.120643, -0.5659155}, // 54799
  {0.056232, 0.120581, -0.5662984}, // 54800
  {0.052411, 0.120576, -0.5667313}, // 54801
  {0.048591, 0.120629, -0.5672751}, // 54802
  {0.044773, 0.120739, -0.5679754}, // 54803
  {0.040958, 0.120908, -0.5688561}, // 54804
  {0.037147, 0.121134, -0.5699194}, // 54805
  {0.033341, 0.121418, -0.5711430}, // 54806
  {0.029541, 0.121759, -0.5724719}, // 54807
  {0.025748, 0.122158, -0.5739260}, // 54808
  {0.021963, 0.122615, -0.5752938}, // 54809
  {0.018187, 0.123129, -0.5764801}, // 54810
  {0.014421, 0.123700, -0.5774336}, // 54811
  {0.010665, 0.124328, -0.5781753}, // 54812
  {0.006921, 0.125013, -0.5788010}, // 54813
  {0.003190, 0.125755, -0.5794499}, // 54814
  {-0.000528, 0.126554, -0.5802526}, // 54815
  {-0.004231, 0.127409, -0.5812807}, // 54816
  {-0.007919, 0.128320, -0.5825247}, // 54817
  {-0.011590, 0.129287, -0.5839062}, // 54818
  {-0.015244, 0.130310, -0.5853138}, // 54819
  {-0.018880, 0.131389, -0.5866413}, // 54820
  {-0.022496, 0.132522, -0.5878120}, // 54821
  {-0.026092, 0.133711, -0.5887843}, // 54822
  {-0.029668, 0.134954, -0.5895486}, // 54823
  {-0.033221, 0.136252, -0.5901207}, // 54824
  {-0.036752, 0.137604, -0.5905379}, // 54825
  {-0.040259, 0.139009, -0.5908556}, // 54826
  {-0.043741, 0.140468, -0.5911419}, // 54827
  {-0.047198, 0.141980, -0.5914685}, // 54828
  {-0.050629, 0.143544, -0.5918997}, // 54829
  {-0.054032, 0.145161, -0.5924836}, // 54830
  {-0.057407, 0.146830, -0.5932450}, // 54831
};

const int s2_mjd_range[2] = { 54618, 54831 };

int
get_iers_bulletin_a( int mjd, float *dut1, float *x, float *y )
{
    if ( (s1_mjd_range[0] <= mjd) && (mjd <= s1_mjd_range[1]) )
    {
        int k = mjd - s1_mjd_range[0];
        *dut1 = s1[k].dut1;
        *x = s1[k].x;
        *y = s1[k].y;
    }
    else if ( (s2_mjd_range[0] <= mjd) && (mjd <= s2_mjd_range[1]) )
    {
        int k = mjd - s2_mjd_range[0];
        *dut1 = s2[k].dut1;
        *x = s2[k].x;
        *y = s2[k].y;
    }
    else
    {
      printf("iers bulletin out of range.  setting to zero.\n");
      *x=0;
      *y=0;
      *dut1=0;
      //assert( 1 == 0 );
      //return -1;
      return 0;
    }
    
    return 0;
}