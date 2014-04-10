function readLoads( path )
    loads_df = readdlm(path, ',')
    #burn the first row and first column
    loads_df = loads_df[2:end, 2:end]
    dts = map(s->strip(s, '"'), loads_df[:, 1])
    vals = loads_df[:, 2:end] * 1e-3
    return dts, vals
end
