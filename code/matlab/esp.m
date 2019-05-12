function [x_esp, y_esp] = esp(img1, img2, options)

    sz = size(img1);
    assert size(img2) == sz;
    d = options.direction;

    [refF2, errormapF2] = get_kSpaceMaps(img1, img2);

    switch options.parametrization
      case 'shell'
        [x, f, y] = esp_shell(refF2, errormapF2);
      case 'cone'
        [x, f, y] = esp_cone(refF2, errormapF2, d);
      case 'zeroConeDistance'
        [x, f, y] = esp_zcd(refF2, errormapF2, d);
    end

    [x_esp, y_esp] = get_espFrac(x, f, y, options.smoothing);

end