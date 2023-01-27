from uproot import concatenate

if __name__=="__main__":
    LAr_position=[-2800.,-2575.,-2400.,-2175.,-2000.,-1775.,-1600.,-1375.,-1200.,-975.,-800.,-575.,-400.,-175.,0.]
    vertex_position=[-299.,-292.,-285.,-278.,-271.,-264.,-216.,-168.,-120.,-72.,-24.,24.,72.,120.,168.,216.,264.,
                    271.,278.,285.,292.,299.]
    effTree=concatenate("/storage/shared/barwu/10thTry/FDEff/FDGeoEff_62877585_990_Eff.root:event_data",
                        ["hadron_selected_eff"], library="np")
    for i_event in range(len(effTree['hadron_selected_eff'])):
        for det_pos in range(len(effTree['hadron_selected_eff'][i_event])):
            for vtx_pos in range(len(effTree['hadron_selected_eff'][i_event][det_pos])):
                print("event ",end="#")
                print(i_event,end=", ")
                print("LAr position",end="=")
                print(LAr_position[det_pos],end=", ")
                print("vertex position",end="=")
                print(vertex_position[vtx_pos],end=", ")
                print("hadron efficiency:",end=" ")
                print(effTree['hadron_selected_eff'][i_event][det_pos][vtx_pos])