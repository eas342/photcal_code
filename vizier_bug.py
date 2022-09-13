# coding: utf-8
from astroquery.vizier import Vizier
Vizier.ROW_LIMIT = 10
tm_res_all = Vizier(catalog='II/246/out').query_region("NGC 2506",radius='0.3deg',get_query_payload=True)
tm_res_all
