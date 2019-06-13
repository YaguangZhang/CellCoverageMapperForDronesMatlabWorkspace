addpath('../')

region = fetchregion([33.701053, 33.813230], [-118.433046, -118.281297], ...
                     'display', true);
elevData = region.readelevation([33.701053, 33.813230], ...
                                [-118.433046, -118.281297], ...
                                'SampleFactor', 10, ...
                                'display', true);
dispelev(elevData, 'mode', 'cartesian');