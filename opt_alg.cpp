#include "opt_alg.h"

solution
MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        while (true) {
            Xopt = rand_mat(N);
            for (int i = 0; i < N; ++i)
                Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
            Xopt.fit_fun(ff, ud1, ud2);
            if (Xopt.y < epsilon) {
                Xopt.flag = 1;
                break;
            }
            if (solution::f_calls > Nmax) {
                Xopt.flag = 0;
                break;
            }
        }
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution MC(...):\n" + ex_info);
    }
}

double *
expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2) {
    try {
        auto *p = new double[2]{0, 0};
        int i = 0;

        solution X0(x0), X1(x0 + d);
        X0.fit_fun(ff, ud1, ud2);
        X1.fit_fun(ff, ud1, ud2);

        if (X0.y == X1.y) {
            p[0] = m2d(X0.x);
            p[1] = m2d(X1.x);
            return p;
        }

        if (X1.y > X0.y) {
            d *= -1;
            X1.x = X0.x + d;
            X1.fit_fun(ff, ud1, ud2);
            if (X1.y >= X0.y) {
                p[0] = m2d(X1.x);
                p[1] = m2d(X0.x - d);
                return p;
            }
        }

        solution X2;
        do {
            X2.x = X0.x + pow(alpha, i) * d;
            X2.fit_fun(ff, ud1, ud2);

            if (solution::f_calls > Nmax) {
                p[0] = m2d(X0.x);
                p[1] = m2d(X2.x);
                return p;
            }

            if (X2.y >= X1.y)
                break;

            X0 = X1;
            X1 = X2;
            i++;
        } while (true);

        if (d > 0) {
            p[0] = m2d(X0.x);
            p[1] = m2d(X2.x);
        } else {
            p[0] = m2d(X2.x);
            p[1] = m2d(X0.x);
        }

        return p;
    }
    catch (string ex_info) {
        throw ("double* expansion(...):\n" + ex_info);
    }
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        unsigned int k = 1;

        vector<double> fibonacci;
        fibonacci.push_back(1);
        fibonacci.push_back(1);

        while (fibonacci[k] < (b - a) / epsilon) {
            fibonacci.push_back(fibonacci[k - 1] + fibonacci[k]);
            k++;
        }

        double c = b - fibonacci[k - 1] / fibonacci[k] * (b - a);
        solution C(c);
        C.fit_fun(ff, ud1, ud2);

        double d = a + b - c;
        solution D(d);
        D.fit_fun(ff, ud1, ud2);

        for (int i = 0; i < k - 3; i++) {
            if (C.y < D.y)
                b = d;
            else
                a = c;

            c = b - fibonacci[k - i - 2] / fibonacci[k - i - 1] * (b - a);
            C = solution(c);
            C.fit_fun(ff, ud1, ud2);

            d = a + b - c;
            D = solution(d);
            D.fit_fun(ff, ud1, ud2);
        }
        Xopt = solution(c);
        Xopt.fit_fun(ff, ud1, ud2);
        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution fib(...):\n" + ex_info);
    }

}

solution
lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1,
    matrix ud2) {
    try {
        solution Xopt;

        solution A(a), B(b);
        double c = (a + b) / 2;
        solution C(c);

        A.fit_fun(ff, ud1, ud2);
        B.fit_fun(ff, ud1, ud2);
        C.fit_fun(ff, ud1, ud2);

        double d_prev = INFINITY;
        int iter = 0;

        while (true) {
            double l = m2d(A.y * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y * (pow(C.x(0), 2) - pow(A.x(0), 2)) +
                           C.y * (pow(A.x(0), 2) - pow(B.x(0), 2)));
            double m = m2d(A.y * (m2d(B.x) - m2d(C.x)) + B.y * (m2d(C.x) - m2d(A.x)) + C.y * (m2d(A.x) - m2d(B.x)));

            if (abs(m) <= epsilon) {
                Xopt = C;
                Xopt.flag = -1;
                return Xopt;
            }

            double d = 0.5 * (l / m);

            if (d < min(m2d(A.x), m2d(B.x)) || d > max(m2d(A.x), m2d(B.x))) {
                Xopt.flag = -2;
                return Xopt;
            }

            solution D(d);
            D.fit_fun(ff, ud1, ud2);

            if (solution::f_calls > Nmax) {
                Xopt = D;
                Xopt.flag = 0;
                return Xopt;
            }

            if (abs(d - d_prev) < gamma) {
                Xopt = D;
                Xopt.flag = 1;
                return Xopt;
            }

            if (d < m2d(C.x)) {
                if (m2d(D.y) < m2d(C.y)) {
                    B = C;
                    C = D;
                } else {
                    A = D;
                }
            } else {
                if (m2d(D.y) < m2d(C.y)) {
                    A = C;
                    C = D;
                } else {
                    B = D;
                }
            }

            double new_interval = abs(m2d(B.x) - m2d(A.x));
            Xopt.ud.add_row(new_interval);

            if (new_interval < epsilon) {
                Xopt = C;
                Xopt.flag = 1;
                return Xopt;
            }

            d_prev = d;
            iter++;
        }
    } catch (...) {
        throw "Error in lag";
    }
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        // Otwieranie pliku do zapisu iteracji
        std::ofstream logFile("HJ_Iterations.csv", std::ios::out);
        logFile << "Iteration;X0;X1;Func_Value\n";

        solution Xopt(x0);
        Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
        solution xB = Xopt;
        int f_calls = 0;

        int iteration = 0;
        while (s > epsilon) {
            solution x_trial = HJ_trial(ff, xB, s, ud1, ud2);

            if (x_trial.y < xB.y) {
                do {
                    xB = x_trial;
                    matrix new_point = xB.x * 2 - Xopt.x;
                    x_trial = solution(new_point);
                    x_trial.y = x_trial.fit_fun(ff, ud1, ud2);
                    f_calls++;

                    // Zapis bie¿¹cej iteracji wewnêtrznej do pliku
                    logFile << iteration << ";" << xB.x(0, 0) << ";" << xB.x(1, 0) << ";" << xB.y << "\n";

                    if (f_calls > Nmax) throw std::runtime_error("Exceeded max function calls");

                } while (x_trial.y < xB.y);

                Xopt = xB;
            } else {
                s *= alpha;
            }

            // Zapis bie¿¹cej iteracji zewnêtrznej do pliku, nawet jeœli x_trial nie zosta³ zaktualizowany
            logFile << iteration << ";" << xB.x(0, 0) << ";" << xB.x(1, 0) << ";" << xB.y << "\n";
            iteration++;

            if (f_calls > Nmax) throw std::runtime_error("Exceeded max function calls");
        }

        logFile.close();  // Zamkniêcie pliku po zakoñczeniu iteracji
        return Xopt;
    } catch (string ex_info) {
        throw ("solution HJ(...):\n" + ex_info);
    }
}


solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2) {
    try {
        solution x_trial = XB;
        int* size = get_size(XB.x);
        int n = size[0];
        delete[] size;

        for (int j = 0; j < n; ++j) {
            matrix step_forward = XB.x;
            step_forward(j, 0) += s;
            solution trial_fwd(step_forward);
            trial_fwd.y = trial_fwd.fit_fun(ff, ud1, ud2);

            if (trial_fwd.y < x_trial.y) {
                x_trial = trial_fwd;
            } else {
                matrix step_backward = XB.x;
                step_backward(j, 0) -= s;
                solution trial_bwd(step_backward);
                trial_bwd.y = trial_bwd.fit_fun(ff, ud1, ud2);

                if (trial_bwd.y < x_trial.y) {
                    x_trial = trial_bwd;
                }
            }
        }

        return x_trial;
    }
    catch (string ex_info) {
        throw ("solution HJ_trial(...):\n" + ex_info);
    }
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        // Otwieranie pliku do zapisu iteracji
        std::ofstream logFile("Rosen_Iterations.csv", std::ios::out);
        logFile << "Iteration;X0;X1;Func_Value\n";

        solution Xopt(x0);
        Xopt.y = Xopt.fit_fun(ff, ud1, ud2);
        matrix x = x0;
        matrix s = s0;
        int iter_count = 0;
        int n = get_len(x0);

        matrix dir = ident_mat(n);

        while (iter_count < Nmax) {
            solution x_trial = Xopt;

            for (int j = 0; j < n; ++j) {
                matrix step = s(j) * dir[j];
                solution forward_trial(x + step);
                forward_trial.y = forward_trial.fit_fun(ff, ud1, ud2);

                solution backward_trial(x - step);
                backward_trial.y = backward_trial.fit_fun(ff, ud1, ud2);

                if (forward_trial.y < Xopt.y) {
                    Xopt = forward_trial;
                    x = Xopt.x;
                    s(j) *= alpha;
                } else if (backward_trial.y < Xopt.y) {
                    Xopt = backward_trial;
                    x = Xopt.x;
                    s(j) *= alpha;
                } else {
                    s(j) *= beta;
                }
            }

            // Zapis bie¿¹cej iteracji do pliku
            logFile << iter_count << ";" << x(0, 0) << ";" << x(1, 0) << ";" << Xopt.y << "\n";
            iter_count++;

            if (norm(s) < epsilon) break;
        }

        logFile.close();  // Zamkniêcie pliku po zakoñczeniu iteracji
        return Xopt;
    } catch (string ex_info) {
        throw ("solution Rosen(...):\n" + ex_info);
    }
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1,
             matrix ud2) {
    try {
        solution Xopt(x0);  // Inicjalizacja rozwi¹zania
        int iter_count = 0;

        while (iter_count < Nmax) {
            iter_count++;

            // Wyznacz funkcjê celu z kar¹
            auto penalized_ff = [&](matrix x, matrix ud1, matrix ud2) {
                matrix f_val = ff(x, ud1, ud2);  // Wartoœæ funkcji celu
                double S = 0.0;  // Kara
                int n = get_len(x);
                for (int i = 0; i < n; i++) {
                    S += std::max(0.0, x(i, 0));  // Oblicz karê (zewnêtrzn¹)
                }
                return f_val + c * S;
            };

            // Optymalizuj funkcjê celu z kar¹ (wykorzystujemy Hooke-Jeeves lub inn¹ metodê optymalizacji)
            Xopt = HJ(penalized_ff, x0, 0.1, 0.5, epsilon, Nmax, ud1, ud2);

            // Skaluj wspó³czynnik kary
            c *= dc;

            // SprawdŸ warunek zakoñczenia
            if (norm(Xopt.x - x0) < epsilon) {
                break;
            }

            // Zaktualizuj punkt startowy
            x0 = Xopt.x;
        }

        return Xopt;
    } catch (string ex_info) {
        throw ("solution pen(...):\n" + ex_info);
    }
}

solution
sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta,
       double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        int n = get_len(x0);  // Wymiar problemu
        std::vector<solution> vertices(n + 1);  // Wierzcho³ki sympleksu
        solution Xopt;  // Rozwi¹zanie optymalne

        // Inicjalizacja sympleksu
        vertices[0] = solution(x0);
        for (int i = 1; i <= n; i++) {
            matrix ei = ident_mat(n)[i - 1];  // Jednostkowy kierunek
            vertices[i] = solution(x0 + s * ei);
        }

        int iter_count = 0;

        while (iter_count < Nmax) {
            iter_count++;

            // Oblicz wartoœci funkcji celu dla wszystkich wierzcho³ków
            for (auto &vertex : vertices) {
                vertex.y = vertex.fit_fun(ff, ud1, ud2);
            }

            // ZnajdŸ indeksy pmin i pmax
            int pmin = 0, pmax = 0;
            for (int i = 1; i <= n; i++) {
                if (vertices[i].y < vertices[pmin].y) pmin = i;
                if (vertices[i].y > vertices[pmax].y) pmax = i;
            }

            // Oblicz œrodek ciê¿koœci (bez pmax)
            matrix p = matrix(n, 1, 0.0);
            for (int i = 0; i <= n; i++) {
                if (i != pmax) p = p + vertices[i].x;
            }
            p = p / n;

            // Odbicie
            matrix podb_x = p + alpha * (p - vertices[pmax].x);
            solution podb(podb_x);
            podb.y = podb.fit_fun(ff, ud1, ud2);

            if (podb.y < vertices[pmin].y) {
                // Ekspansja
                matrix pe_x = p + gamma * (podb.x - p);
                solution pe(pe_x);
                pe.y = pe.fit_fun(ff, ud1, ud2);

                if (pe.y < podb.y) {
                    vertices[pmax] = pe;
                } else {
                    vertices[pmax] = podb;
                }
            } else if (podb.y < vertices[pmax].y) {
                vertices[pmax] = podb;
            } else {
                // Zawê¿enie
                matrix pz_x = p + beta * (vertices[pmax].x - p);
                solution pz(pz_x);
                pz.y = pz.fit_fun(ff, ud1, ud2);

                if (pz.y < vertices[pmax].y) {
                    vertices[pmax] = pz;
                } else {
                    // Redukcja
                    for (int i = 0; i <= n; i++) {
                        if (i != pmin) {
                            vertices[i].x = delta * (vertices[i].x + vertices[pmin].x);
                            vertices[i].y = vertices[i].fit_fun(ff, ud1, ud2);
                        }
                    }
                }
            }

            // Warunek zakoñczenia
            double max_dist = 0.0;
            for (int i = 0; i <= n; i++) {
                max_dist = std::max(max_dist, norm(vertices[pmin].x - vertices[i].x));
            }
            if (max_dist < epsilon) break;
        }

        Xopt = vertices[pmin];  // Najlepszy wierzcho³ek
        return Xopt;
    } catch (string ex_info) {
        throw ("solution sym_NM(...):\n" + ex_info);
    }
}

solution
SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution SD(...):\n" + ex_info);
    }
}

solution
CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution CG(...):\n" + ex_info);
    }
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
                matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1,
                matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Newton(...):\n" + ex_info);
    }
}

solution
golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution golden(...):\n" + ex_info);
    }
}

solution Powell(matrix(*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution Powell(...):\n" + ex_info);
    }
}

solution
EA(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon,
   int Nmax, matrix ud1, matrix ud2) {
    try {
        solution Xopt;
        //Tu wpisz kod funkcji

        return Xopt;
    }
    catch (string ex_info) {
        throw ("solution EA(...):\n" + ex_info);
    }
}
