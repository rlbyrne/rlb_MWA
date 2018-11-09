pro catalog_manipulation

  restore, '/Users/rubybyrne/FHD/catalog_data/master_sgal_cat.sav', /relaxed
  fluxes = []
  for i=0,n_elements(catalog)-1 do fluxes=[fluxes, catalog[i].flux.i]
  catalog_sorted = catalog[reverse(sort(fluxes))]
  catalog_4k = catalog_sorted[0:3999]  ; choose the brightest 4,000 sources to eliminate spatial dependence
  random_nums = sort(randomu(seed, n_elements(catalog_4k)))  ; sorting random numbers gives random unique indices
  catalog_2p5k = catalog_sorted[random_nums[0:2499]]
  catalog_2k = catalog_sorted[random_nums[0:1999]]
  catalog_1p5k = catalog_sorted[random_nums[0:1499]]
  neg_sources = catalog_sorted[random_nums[2000:2499]]
  for i=0,n_elements(neg_sources)-1 do neg_sources[i].flux.I = -neg_sources[i].flux.I
  catalog_2p5k_neg_sources = [catalog_sorted[random_nums[0:1999]], neg_sources]
  save, catalog_2p5k, filename='/Users/rubybyrne/incomplete_sky_models_for_testing/test_catalog_2p5k_sources.sav'
  save, catalog_2k, filename='/Users/rubybyrne/incomplete_sky_models_for_testing/test_catalog_2k_sources.sav'
  save, catalog_1p5k, filename='/Users/rubybyrne/incomplete_sky_models_for_testing/test_catalog_1p5k_sources.sav'
  save, catalog_2p5k_neg_sources, filename='/Users/rubybyrne/incomplete_sky_models_for_testing/test_catalog_2p5k_neg_sources.sav'

end