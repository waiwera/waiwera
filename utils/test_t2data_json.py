import unittest
from t2data_json import *

class t2data_jsonTestCase(unittest.TestCase):

    def test_generators(self):
        """Test conversion of generators to sources"""

        geo = mulgrid().rectangular([100.]*10, [150.]*10, [20.]*10)
        dat = t2data_export_json()
        eosname = 'wce'

        # Constant mass production:
        name, blkindex, q = 'prd01', 5, -12.3
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'MASS', gx = q)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['rate'], q)

        # Constant mass injection:
        name, blkindex, q, h = 'inj01', 12, 5.9, 1800.
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'MASS', gx = q, ex = h)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['rate'], q)
        self.assertEqual(json['source'][-1]['enthalpy'], h)
        
        # Constant heat injection:
        name, blkindex, q = 'hea01', 10, 2400.
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'HEAT', gx = q)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['rate'], q)
        self.assertEqual(json['source'][-1]['component'], 3)

        # DELV:
        name, blkindex, PI, Pb = 'del01', 10, 1.e-12, 1.5e5
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELV', gx = PI, ex = Pb)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['deliverability']['productivity'], PI)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)

        # DELG:
        name, blkindex, PI, Pb, limit = 'del02', 11, 1.e-12, 1.5e5, None
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELG', gx = PI, ex = Pb, hg = limit)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['deliverability']['productivity'], PI)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)
        self.assertEqual(json['source'][-1]['direction'], 'production')

        # DELG with steam limiter:
        name, blkindex, PI, Pb, limit, Ps = 'del03', 13, 1.e-11, 2.5e5, 16., None
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELG', gx = PI, ex = Pb, hg = limit, fg = Ps)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['deliverability']['productivity'], PI)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)
        self.assertEqual(json['source'][-1]['direction'], 'production')
        self.assertEqual(json['source'][-1]['limiter']['type'], 'steam')
        self.assertEqual(json['source'][-1]['limiter']['limit'], limit)

        # DELG with steam limiter and separator pressure:
        name, blkindex, PI, Pb, limit, Ps = 'del04', 14, 1.e-11, 2.5e5, 16., 2.e5
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELG', gx = PI, ex = Pb, hg = limit, fg = Ps)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['deliverability']['productivity'], PI)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)
        self.assertEqual(json['source'][-1]['direction'], 'production')
        self.assertEqual(json['source'][-1]['limiter']['type'], 'steam')
        self.assertEqual(json['source'][-1]['limiter']['limit'], limit)
        self.assertEqual(json['source'][-1]['limiter']['separator_pressure'], Ps)

        # DELG with initial rate specified (to calculate PI):
        name, blkindex, PI, Pb, rate, Ps = 'del05', 14, 1.e-11, 2.5e5, -8., 2.e5
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELG', gx = PI, ex = Pb, hg = rate, fg = Ps)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['rate'], rate)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)
        self.assertEqual(json['source'][-1]['direction'], 'production')
        self.assertFalse('productivity' in json['source'][-1]['deliverability'])

        # DELS with steam limiter and separator pressure:
        name, blkindex, PI, Pb, limit, Ps = 'del06', 14, 1.e-11, 2.5e5, 16., 2.e5
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELS', gx = PI, ex = Pb, hg = limit, fg = Ps)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['deliverability']['productivity'], PI)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)
        self.assertEqual(json['source'][-1]['direction'], 'production')
        self.assertEqual(json['source'][-1]['limiter']['type'], 'steam')
        self.assertEqual(json['source'][-1]['limiter']['limit'], limit)
        self.assertEqual(json['source'][-1]['limiter']['separator_pressure'], Ps)

        # DELW with steam limiter and separator pressure:
        name, blkindex, PI, Pb, limit, Ps = 'del07', 15, 1.e-11, 2.5e5, 16., 2.e5
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELW', gx = PI, ex = Pb, hg = limit, fg = Ps)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['deliverability']['productivity'], PI)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)
        self.assertEqual(json['source'][-1]['direction'], 'production')
        self.assertEqual(json['source'][-1]['limiter']['type'], 'water')
        self.assertEqual(json['source'][-1]['limiter']['limit'], limit)
        self.assertEqual(json['source'][-1]['limiter']['separator_pressure'], Ps)

        # DELT:
        name, blkindex, PI, Pb, limit, Ps = 'del08', 16, 1.e-11, 2.5e5, 16., 2.e5
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'DELT', gx = PI, ex = Pb, hg = limit, fg = Ps)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['deliverability']['productivity'], PI)
        self.assertEqual(json['source'][-1]['deliverability']['pressure'], Pb)
        self.assertEqual(json['source'][-1]['direction'], 'production')
        self.assertEqual(json['source'][-1]['limiter']['type'], 'total')
        self.assertEqual(json['source'][-1]['limiter']['limit'], limit)
        self.assertFalse('separator_pressure' in json['source'][-1]['limiter'])

        # RECH:
        name, blkindex, rate, h, hg, dirn = 'del09', 17, 13.5, 2.5e5, None, None
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'RECH', gx = rate, ex = h, hg = hg, fg = dirn)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['rate'], rate)
        self.assertEqual(json['source'][-1]['enthalpy'], h)

        # RECH with pressure reference:
        name, blkindex, coef, h, Pr, dirn = 'del10', 18, 1.e-12, 1500., 2.5e5, None
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'RECH', gx = coef, ex = h, hg = Pr, fg = dirn)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['enthalpy'], h)
        self.assertEqual(json['source'][-1]['direction'], 'both')
        self.assertEqual(json['source'][-1]['recharge']['pressure'], Pr)
        self.assertEqual(json['source'][-1]['recharge']['coefficient'], coef)

        # RECH with pressure reference and direction:
        name, blkindex, coef, h, Pr, dirn = 'del11', 19, 1.e-12, 1500., 2.5e5, -1.0
        gen = t2generator(name, geo.block_name_list[blkindex],
                          type = 'RECH', gx = coef, ex = h, hg = Pr, fg = dirn)
        dat.add_generator(gen)
        json = dat.generators_json(geo, eosname = eosname)
        self.assertEqual(json['source'][-1]['name'], name)
        self.assertEqual(json['source'][-1]['cell'], blkindex)
        self.assertEqual(json['source'][-1]['enthalpy'], h)
        self.assertEqual(json['source'][-1]['direction'], 'out')
        self.assertEqual(json['source'][-1]['recharge']['pressure'], Pr)
        self.assertEqual(json['source'][-1]['recharge']['coefficient'], coef)
        
if __name__ == '__main__':

    suite = unittest.TestLoader().loadTestsFromTestCase(t2data_jsonTestCase)
    unittest.TextTestRunner(verbosity = 1).run(suite)
        
