
template<unsigned int N>
void vtx_pack(TkObj track[N], MP7DataWord data[]){
    for(unsigned int i = 0; i < N; i++){
        data[i] = 0;
        data[i](7, 0) = VxfZ0(track[i]);
        data[i](23, 8) = track[i].hwPt;
        if(track[i].hwPt > 0){
            data[i][24] = 1;
        }else{
            data[i][24] = 0;
        }
    }
}

template<unsigned int N>
void vxf_pattern_file(MP7PatternSerializer & vtxPatterns, TkObj track[N]){
    //MP7PatternSerializer vtxPatterns("vtx_patterns.txt", 1, 0);
    // For simplicity, initialise MP7_NCHANN channels
    // But the channels unused by vertexing finder are always 0
    MP7DataWord frame[MP7_NCHANN];
    for(int ic = VTX_NCHANN; ic < MP7_NCHANN; ic++){
        frame[ic] = 0;
    }

    TkObj trackFrame[VTX_NCHANN];
    int it = 0;
    while(it < N){
        for(int ic = 0; ic < VTX_NCHANN; ic++){
            if(it < N){
                trackFrame[ic] = track[it];
            }else{
                clear(trackFrame[ic]);
            }
            it++;
        }
        vtx_pack<VTX_NCHANN>(trackFrame, frame);
        vtxPatterns(frame);
    }

    for(int ic = VTX_NCHANN; ic < MP7_NCHANN; ic++){
        frame[ic] = 0;
    }
    // Inter-event gap
    for(it = 0; it < 40; it++){
        vtxPatterns.print(0, frame, true);
    }
}
