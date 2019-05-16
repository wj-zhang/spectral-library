import unittest
import fragments
from modifications import modifications_env


class FragmentsTestCase(unittest.TestCase):
    """Tests for 'fragments.py'."""

    def setUp(self):
        self.peptide = 'RPLEDGDQPDAK'
        self.masses = [95.02481842041016, 97.02877807617188, 98.06063842773438, 98.98450469970703, 101.07130432128906, 102.05528259277344, 104.05342864990234, 110.07156372070312, 111.04438018798828, 112.05085754394531, 112.08729553222656, 113.03451538085938, 113.05418395996094, 113.58711242675781, 115.08719635009766, 120.08098602294922, 127.05017852783203, 129.06607055664062, 129.10240173339844, 130.05001831054688, 130.08644104003906, 133.06085205078125, 136.06179809570312, 136.07620239257812, 137.06541442871094, 138.03746032714844, 139.15802001953125, 145.0607147216797, 147.1127471923828, 147.1234130859375, 149.02337646484375, 149.04501342773438, 150.0447540283203, 151.04153442382812, 155.08151245117188, 157.10787963867188, 159.0766143798828, 167.03411865234375, 167.05519104003906, 167.08106994628906, 168.0557098388672, 169.05279541015625, 172.07180786132812, 173.0558624267578, 175.1192169189453, 179.89044189453125, 181.06088256835938, 181.0973358154297, 183.14891052246094, 185.09222412109375, 187.0715789794922, 194.12933349609375, 195.11268615722656, 198.08775329589844, 199.0712432861328, 200.1393585205078, 201.12391662597656, 208.95318603515625, 209.09249877929688, 209.1392822265625, 210.9508514404297, 211.14341735839844, 213.08692932128906, 215.13934326171875, 215.6182403564453, 216.09861755371094, 216.12034606933594, 218.14991760253906, 219.15304565429688, 223.06362915039062, 224.06643676757812, 225.04327392578125, 226.0435333251953, 226.0826873779297, 227.0415496826172, 227.0662078857422, 232.0446014404297, 237.13548278808594, 243.13357543945312, 244.09222412109375, 245.07620239257812, 248.64724731445312, 250.0943145751953, 254.16107177734375, 255.16456604003906, 261.15753173828125, 263.1467590332031, 264.2062683105469, 267.1621398925781, 267.66668701171875, 270.0720520019531, 270.64154052734375, 283.103515625, 284.0890808105469, 284.12408447265625, 288.0832214355469, 290.184814453125, 297.15802001953125, 298.637939453125, 301.11407470703125, 302.1337890625, 305.107421875, 308.1943359375, 309.1173400878906, 315.16583251953125, 318.9222106933594, 320.9189147949219, 321.2391052246094, 323.1322021484375, 324.1189880371094, 325.6670227050781, 326.9854431152344, 333.17584228515625, 334.1813049316406, 334.6729736328125, 335.0975646972656, 339.2502136230469, 340.2558288574219, 341.0191650390625, 342.0196228027344, 343.0165100097656, 344.0168762207031, 344.9939880371094, 345.99603271484375, 348.0698547363281, 349.2341003417969, 350.23681640625, 353.1097412109375, 358.16265869140625, 359.0281982421875, 360.028564453125, 361.0263366699219, 362.0259094238281, 367.1728820800781, 367.2451171875, 368.24884033203125, 374.174560546875, 383.18145751953125, 392.1893005371094, 393.26214599609375, 394.2104187011719, 399.11639404296875, 412.21868896484375, 429.085205078125, 430.2296142578125, 430.8869934082031, 431.0889587402344, 431.2328186035156, 439.73577880859375, 439.9068908691406, 440.0718078613281, 440.2412414550781, 447.0361022949219, 448.234375, 451.26416015625, 468.2920837402344, 469.2941589355469, 476.11676025390625, 479.2586669921875, 487.2565002441406, 487.6568298339844, 496.2887268066406, 497.2889709472656, 513.4737548828125, 518.4730224609375, 518.674072265625, 518.8762817382812, 524.2772216796875, 539.3308715820312, 540.27490234375, 541.2619018554688, 542.2670288085938, 554.2821044921875, 555.03662109375, 558.2849731445312, 569.0481567382812, 585.1800537109375, 590.9908447265625, 593.3049926757812, 594.292236328125, 611.3156127929688, 612.3198852539062, 636.3226928710938, 636.8225708007812, 650.3232421875, 651.3158569335938, 668.3362426757812, 669.3391723632812, 673.3111572265625, 682.3434448242188, 694.35009765625, 711.378173828125, 730.3358154296875, 765.3524169921875, 766.3485717773438, 767.3388061523438, 782.3684692382812, 783.3638305664062, 784.3653564453125, 800.3924560546875, 801.3890380859375, 845.3566284179688, 911.4097900390625, 1007.2291259765625, 1100.7218017578125, 1349.7861328125]
        self.intensities = [3643.9434, 20864.805, 3150.0356, 49744.945, 64602.81, 68801.48, 4412.7334, 22118.564, 44962.734, 192638.7, 10895.325, 5446.2583, 4734.6567, 10111.674, 3264.7612, 5891.307, 20858.12, 14974.961, 94549.84, 17911.385, 32727.717, 20196.172, 111258.0, 5422.9023, 4256.448, 2546.054, 2974.4187, 17152.312, 74371.83, 4464.506, 321687.47, 28999.395, 5476.986, 6483.0146, 3086.2214, 19398.926, 3580.1982, 54340.42, 42226.04, 3102.2964, 5613.659, 18115.562, 4145.3965, 42507.43, 15668.241, 3489.3564, 4768.943, 4564.3545, 4512.3486, 12217.58, 3883.8137, 3437.403, 16792.494, 13498.8545, 6193.7446, 21615.99, 13209.15, 52366.633, 14008.867, 2964.1921, 17658.348, 4895.6943, 57714.95, 12451.928, 43622.902, 5593.5596, 4461.6753, 197499.88, 11735.32, 16908.19, 5390.81, 14485.062, 3528.539, 15523.082, 9859.136, 19253.56, 4405.147, 17460.494, 6532.5415, 11378.965, 4630.579, 3294.2068, 4828.175, 118236.67, 10962.832, 3812.1045, 2898.4758, 3648.09, 3949.5027, 3070.5586, 5657.093, 11919.075, 4284.904, 11511.396, 38335.984, 11329.403, 3536.4849, 3584.7053, 5374.4795, 3200.5896, 5951.0146, 2722.5154, 4297.3794, 4451.275, 12504.423, 16684.531, 45292.8, 14524.298, 5983.292, 20768.275, 4706.1616, 6529.6226, 47741.18, 4383.081, 4588.278, 10806.4375, 25399.184, 4827.058, 11966.245, 5494.644, 34184.78, 6274.234, 12496.1455, 3968.2432, 17376.87, 36977.203, 4850.3257, 6052.38, 4219.186, 62119.594, 64641.03, 308298.5, 66173.66, 3974.4937, 150304.11, 29224.795, 3645.2715, 4132.364, 6051.4253, 3613.9182, 3438.477, 5442.889, 17853.494, 3700.4773, 207404.89, 15246.081, 20952.443, 33197.113, 21802.13, 21673.043, 20568.994, 5723.744, 16438.611, 11495.006, 6269.2627, 32783.254, 9968.02, 5169.2686, 4711.9976, 3242.08, 3184.9248, 82282.8, 17836.838, 3504.1282, 14036.317, 17416.494, 14268.485, 4292.7344, 5190.252, 24080.988, 93967.56, 20439.38, 3039.4353, 3232.2544, 4688.1733, 3830.5266, 3471.5627, 3795.3484, 13775.937, 10869.184, 70458.85, 18349.666, 5349.588, 4932.8584, 17038.887, 6198.451, 66830.086, 20706.625, 16651.305, 3594.3813, 4522.6025, 5571.427, 13063.42, 37596.023, 17205.088, 3955.1692, 23781.062, 81719.37, 28056.504, 29365.574, 5367.863, 12153.936, 5864.6543, 3277.7224, 3845.1772, 4128.8716]
        self.error = 0.02

        mod_digits = [4, 35, 7]
        self.mod_env = modifications_env(mod_digits)

        self.modified_peptide = 'DhNGNGTYSeCSYVPR'
        self.modified_masses = [112.08671569824219, 113.21510314941406, 115.08621215820312, 116.07046508789062, 127.04966735839844, 129.10235595703125, 129.1094970703125, 130.08592224121094, 133.0602264404297, 136.07577514648438, 140.08126831054688, 143.04490661621094, 144.0767059326172, 145.0600128173828, 147.1127471923828, 155.0445556640625, 158.09205627441406, 159.07754516601562, 167.04443359375, 169.0974578857422, 169.13304138183594, 171.11256408691406, 172.07107543945312, 173.05618286132812, 175.11880493164062, 178.05303955078125, 181.09706115722656, 185.05577087402344, 190.08238220214844, 191.1025390625, 193.09649658203125, 197.12672424316406, 199.10787963867188, 200.11074829101562, 201.08798217773438, 201.1229705810547, 202.0883331298828, 202.66578674316406, 203.065185546875, 210.02561950683594, 210.08895874023438, 212.1393280029297, 213.0499725341797, 213.12384033203125, 216.13177490234375, 221.09141540527344, 222.094970703125, 223.10699462890625, 224.06671142578125, 225.05001831054688, 226.11834716796875, 226.13674926757812, 227.10191345214844, 227.12191772460938, 228.0630645751953, 228.09814453125, 229.09365844726562, 229.11639404296875, 230.0769500732422, 231.04217529296875, 231.06077575683594, 232.13864135742188, 234.14405822753906, 235.1459197998047, 237.12405395507812, 238.0824737548828, 242.07730102539062, 243.1441650390625, 244.13101196289062, 248.070068359375, 250.1204071044922, 251.10165405273438, 252.0615234375, 252.1324920654297, 254.16082763671875, 255.10906982421875, 255.1448974609375, 256.09130859375, 266.1486511230469, 267.10858154296875, 268.1285400390625, 269.08721923828125, 270.07159423828125, 270.14410400390625, 271.0540771484375, 271.1470947265625, 272.1717224121094, 273.1748962402344, 283.05572509765625, 285.119140625, 287.0977478027344, 288.08306884765625, 288.1552734375, 289.15863037109375, 295.1034851074219, 304.1290588378906, 309.0828857421875, 313.1133117675781, 325.1138000488281, 326.1092834472656, 336.1234130859375, 339.09271240234375, 344.1213073730469, 347.23101806640625, 348.2323303222656, 357.10406494140625, 364.12396240234375, 367.0882873535156, 371.2397155761719, 372.2433776855469, 378.10345458984375, 382.136962890625, 385.100341796875, 392.119873046875, 393.1244201660156, 394.1251525878906, 396.1167297363281, 402.12750244140625, 407.158447265625, 410.1318664550781, 411.1343078613281, 413.18695068359375, 422.2406311035156, 424.1112976074219, 427.1570129394531, 428.1463623046875, 441.1354675292969, 442.12322998046875, 449.20587158203125, 450.20745849609375, 452.1626281738281, 459.1477355957031, 462.2554016113281, 463.2618103027344, 469.24169921875, 479.15289306640625, 480.15570068359375, 484.2800598144531, 488.1731872558594, 489.1332092285156, 490.14031982421875, 496.1830139160156, 497.1598815917969, 498.166259765625, 500.9794921875, 505.2034912109375, 506.165283203125, 507.1484375, 508.1515808105469, 515.1725463867188, 524.1780395507812, 525.1588745117188, 526.1595458984375, 534.3037109375, 535.3048095703125, 537.2465209960938, 542.1856079101562, 543.1743774414062, 544.1715698242188, 554.2726440429688, 557.7510375976562, 558.2551879882812, 560.1948852539062, 563.286865234375, 602.8021240234375, 620.3334350585938, 621.3359375, 622.3388671875, 639.303466796875, 640.2881469726562, 642.2179565429688, 643.0076293945312, 643.2005615234375, 653.1984252929688, 660.2299194335938, 661.2254638671875, 670.2125854492188, 671.796142578125, 677.2498168945312, 678.237548828125, 687.2407836914062, 688.2244873046875, 689.2194213867188, 702.3380126953125, 705.2481689453125, 706.240478515625, 747.2727661132812, 771.4232788085938, 775.2634887695312, 781.3671264648438, 782.3696899414062, 783.373779296875, 792.2772827148438, 796.1099243164062, 796.1648559570312, 796.427978515625, 797.1073608398438, 817.4256591796875, 818.429931640625, 835.4267578125, 850.3887329101562, 868.3973388671875, 869.4026489257812, 870.396484375, 881.4255981445312, 882.4230346679688, 883.4396362304688, 936.4530639648438, 1013.4341430664062, 1031.4610595703125, 1032.4639892578125, 1033.465087890625, 1043.5166015625, 1104.5782470703125, 1106.5426025390625, 1107.5255126953125, 1114.483154296875, 1122.5677490234375, 1128.8607177734375, 1132.50732421875, 1133.5107421875, 1134.514892578125, 1171.52490234375, 1172.50341796875, 1189.531005859375, 1190.533447265625, 1191.534912109375, 1239.664306640625, 1269.5430908203125, 1275.6536865234375, 1302.6636962890625, 1303.6385498046875, 1304.61572265625, 1305.5804443359375, 1326.553955078125, 1343.5751953125, 1344.5712890625, 1361.603759765625, 1403.6595458984375]
        self.modified_intensities = [2430.0405, 1770.075, 7745.551, 2798.1653, 13409.105, 24407.605, 2257.853, 2356.6777, 2282.8477, 28368.535, 2525.5889, 2936.5515, 2414.3254, 2722.48, 4673.132, 3584.0437, 3927.35, 2204.2, 2130.4717, 2475.4607, 16214.956, 18589.775, 20829.28, 3998.3662, 65823.1, 12804.197, 5200.766, 39526.266, 3321.1287, 21321.033, 5802.7935, 3551.9944, 45993.47, 7692.4907, 3151.7092, 2401.9644, 2519.0725, 2167.4365, 3723.6487, 4556.2544, 2678.388, 3280.718, 5542.0522, 3626.1025, 2241.0054, 4017.8215, 2752.548, 3640.3987, 2618.2256, 4020.3462, 32781.887, 3510.2236, 5868.319, 2540.7083, 2423.571, 3590.2324, 4589.0557, 3193.9233, 4732.463, 2243.5615, 25601.926, 2489.0295, 7291.3677, 4791.24, 4411.024, 7371.815, 16611.438, 3682.077, 2920.4226, 4134.6206, 2743.8892, 6408.3984, 4474.354, 2989.0884, 4875.5967, 11925.606, 55896.95, 4893.8896, 3924.717, 5169.582, 4379.323, 4551.13, 19298.346, 12423.599, 3018.712, 3150.5024, 319749.75, 49103.2, 3580.0361, 5818.437, 12980.038, 13003.193, 138556.6, 38845.598, 16012.963, 3782.7312, 3486.355, 6446.3813, 5051.3057, 5208.8228, 3250.4414, 27577.31, 6877.494, 6173.7573, 3431.069, 16174.403, 3609.6934, 22412.11, 75194.836, 15400.012, 6291.442, 3730.8105, 18786.717, 15995.643, 5914.0024, 3456.9456, 7090.955, 18474.268, 3544.3887, 15729.805, 14815.937, 4468.1416, 5696.9707, 12144.9, 4357.7515, 5847.3467, 5948.429, 21672.559, 16831.818, 15737.978, 5145.566, 6310.2485, 4954.161, 7037.325, 7706.3555, 4656.0654, 14039.498, 7572.2915, 3869.4797, 3499.3208, 3831.4993, 4863.8745, 21521.537, 16686.598, 3185.4336, 3886.7607, 3933.5325, 27848.564, 5114.6133, 7613.648, 14535.127, 43919.188, 14232.541, 43022.246, 16389.756, 5150.7236, 24830.38, 18884.668, 3924.4353, 5802.782, 4592.4775, 3312.0188, 5751.319, 4110.282, 4079.662, 6787.3057, 83114.44, 33845.516, 4928.3027, 4569.36, 4625.0054, 3468.0894, 6144.2905, 3727.0664, 12943.842, 12492.389, 3844.169, 4569.8354, 4478.023, 5149.7334, 5698.424, 18074.586, 4265.2266, 4902.0107, 7064.7866, 4723.1543, 4131.738, 4160.512, 6139.8394, 101358.336, 43406.027, 5614.665, 4778.9185, 29521.166, 22276.904, 122951.266, 15312.05, 14822.603, 11125.56, 4899.48, 4818.2446, 259636.8, 148184.6, 36830.39, 17439.512, 20507.502, 5869.3975, 7064.5, 9248.184, 234667.98, 157583.86, 40347.863, 5042.386, 22425.03, 6320.5293, 22180.63, 6819.1484, 6952.0015, 5110.4355, 64472.562, 54761.203, 12115.191, 14348.579, 4590.6274, 95469.45, 87783.84, 15688.616, 25051.895, 4358.65, 7050.2544, 4500.0176, 39779.574, 32015.85, 4942.9976, 5257.7007, 6242.7197, 13000.709, 6020.975, 4424.0005]

    def test_ion_masses(self):
        b_ion_masses, y_ion_masses = fragments.ion_masses(self.peptide)
        self.assertAlmostEqual(b_ion_masses[0], 157.1083, places=2)
        self.assertAlmostEqual(y_ion_masses[-1], 1184.5429, places=2)

        b_ion_masses, y_ion_masses = fragments.ion_masses(self.modified_peptide, aa_mass=self.mod_env['aa_mass'])
        self.assertAlmostEqual(b_ion_masses[0], 116.0342, places=2)
        self.assertAlmostEqual(y_ion_masses[0], 175.1189, places=2)

    def test_ion_annotation(self):
        ion_flags = fragments.ion_annotation(self.masses, self.intensities, self.peptide, error=self.error)

        ion_types = ('b', 'y')
        ion_nums = (1, 2, 3, 4, 5, 6, 7, 8)

        for ion_type in ion_types:
            for ion_num in ion_nums:
                ion_flag = ion_type + str(ion_num)
                self.assertIn(ion_flag, ion_flags)

        ion_flags = fragments.ion_annotation(self.modified_masses, self.modified_intensities, self.modified_peptide,
                                             aa_mass=self.mod_env['aa_mass'], error=self.error)

        ion_flag_keys = ['b2', 'b3', 'b4', 'b5', 'b6', 'y1', 'y2', 'y3', 'y4', 'y5', 'y6', 'y7', 'y8', 'y9', 'y10']
        for ion_flag_key in ion_flag_keys:
            self.assertIn(ion_flag_key, ion_flags)


if __name__ == '__main__':
    unittest.main()
