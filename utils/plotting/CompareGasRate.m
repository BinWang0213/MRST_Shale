function CompareGasRate(ws,ws_r,report)
plotWellSols({ws, ws_r}, report.ReservoirTime,'field','qWs', 'datasetnames', {'EDFM', 'Fully resolved solution'})
end