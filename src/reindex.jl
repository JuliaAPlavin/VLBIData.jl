function uv_reindex(uvtbl)
    if !isempty(uvtbl) && hasoptic(first(uvtbl), antennas)
        antname_to_id = @p uvtbl flatmap(antennas, _2.name) unique StructArray((__, 1:length(__))) Dict
        uvtbl = @modify(antennas(uvtbl[∗])[∗]) do a
            @set a.id = antname_to_id[a.name]
        end
    end
    # XXX: non-composable!
    if !isempty(uvtbl) && hasproperty(first(uvtbl), :freq_spec) && hasproperty(first(uvtbl).freq_spec, :ix)
        freq_to_id = @p uvtbl map(_.freq_spec) unique StructArray((__, 1:length(__))) Dict
        uvtbl = @modify(uvtbl[∗].freq_spec) do f
            @set f.ix = freq_to_id[f]
        end
    end
    return uvtbl
end
