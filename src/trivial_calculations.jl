@stable uvshift(center, uvtbl) = @p uvtbl mapset(value=_.value / visibility(Point(flux=1, coords=center), _.spec.uv))

@stable add_conjvis(uvtbl) = vcat(uvtbl, map(conjvis, uvtbl))
