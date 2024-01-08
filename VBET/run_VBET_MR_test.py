import classVBET


class RunVBET:
    def __init__(self):
        self.params = {
            'network': 'data-raw/flowlines_with_cumarea_clip2.shp',
            'dem': 'data-raw/dem_nhdplushr_yuba_meters_v2.tif',
            'out': 'output/yuba_output_v3.shp',
            'scratch': 'scratch',
            'lg_da': 250,
            'med_da': 25,
            'lg_slope': 3,
            'med_slope': 4,
            'sm_slope': 5,
            'lg_buf': 1000, #500,
            'med_buf': 1000,#200,
            'sm_buf': 1000,#80,
            'min_buf': 10,
            'dr_area': None,
            'da_field': 'tt_d_s_', #'TotDASqKm',
            'lg_depth': 3,
            'med_depth': 2,
            'sm_depth': 1.5
            }

    def run(self):
        vb = classVBET.VBET(**self.params)
        if self.params['da_field'] is None:
            vb.add_da()
        vb.valley_bottom()

vbrun = RunVBET()
vbrun.run()

