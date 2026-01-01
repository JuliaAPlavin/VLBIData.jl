uvshift(uvtbl, center) = @p uvtbl mapset(value=_.value / visibility(Point(flux=1, coords=center), _.spec.uv))

# uvrotate(uvtbl, angle) = @p uvtbl mapset(value=_.value * exp(im * angle))

add_conjvis(uvtbl) = vcat(uvtbl, map(conjvis, uvtbl))
