if (auto fhandler = fopen(argv[1], "r")) {
                        char s[50] = "\0";
                        int comma_indicies[4], stringpos = 0;
                        comma_indicies[4] = -1;

                        if (!fgets(s, 50, fhandler))
                                return 1;


                        for (int com_i = 0; com_i < 4; com_i++) {
                                bool found = false;
                                while(!found and (stringpos < 50)) {
                                        if (s[stringpos] == ',') {
                                                found = true;
                                                comma_indicies[com_i] = stringpos;
                                        }
                                        if ((s[stringpos] == '\0') or (s[stringpos] == EOF))
                                                return 1;
                                        stringpos++;
                                }
                        }

                        if ((comma_indicies[1] - comma_indicies[0] == 1) or
                            (comma_indicies[2] - comma_indicies[1] == 1) or
                            (comma_indicies[3] - comma_indicies[2] == 1))
                                return 1;

                        char buf_nx[10];
                        char buf_ny[10];
                        char buf_k1[10];
                        char buf k2[10];

       
                        printf("LALA %d %d %d %d",
                                        comma_indicies[0],
                                        comma_indicies[1],
                                        comma_indicies[2],
                                        comma_indicies[3]);

//////////////////
vector<bool> initVec(Nx, false);
        vector<vector<bool>> graph(Ny, initVec);
        for (auto j = graph.begin(); j != graph.end(); ++j) {
        for (auto i = j.begin(); i != j.end(); ++i)
                cout << *i << ' ';
        cout << "\n";
        }

