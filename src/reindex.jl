function uv_reindex(uvtbl)
    # XXX: non-composable!
    if !isempty(uvtbl) && hasproperty(first(uvtbl), :freq_spec) && hasproperty(first(uvtbl).freq_spec, :ix)
        freq_to_id = @p uvtbl map(_.freq_spec) unique StructArray((__, 1:length(__))) Dict
        uvtbl = @modify(uvtbl[∗].freq_spec) do f
            @set f.ix = freq_to_id[f]
        end
    end
    return uvtbl
end
