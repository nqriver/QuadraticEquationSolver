#include <iostream>
#include <sstream>
#include <complex>
#include <optional>


template <typename T>
std::ostream &operator<<(std::ostream& out, const std::complex<T> complexNumber) {
    out << complexNumber.real() << " + " << complexNumber.imag() << " * i";
    return out;
}

struct Equation {
private:
    double a, b, c;
    double delta;
public:

    Equation(double a, double b, double c) : a{a}, b{b}, c{c}, delta{ b*b-4*a*c }
    {
    }

    std::string put() const {
        std::ostringstream equation{};
        if (a != 0) equation << a << "x^2";
        if (b != 0) equation << " + " << b << "x";
        if (c != 0) equation << " + " << c << "\n";
        return equation.str();
    }

    void setCoefficients(double a, double b, double c) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->delta = b*b - 4*a*c;
    }

    double getDelta() const {
        return delta;
    }

    double getA() const {
        return a;
    }

    double getB() const {
        return b;
    }

    double getC() const {
        return c;
    }
};


std::ostream &operator<<(std::ostream& out, const Equation& equation) {
    out << equation.put();
    return out;
}

struct EquationSolver {
    using real_solution_t = std::pair<double, double>;
    using complex_solution_t = std::pair<std::complex<double>, std::complex<double>>;

private:
    const Equation &equation;

public:
    explicit EquationSolver(Equation &equation) : equation{equation}
    {
    }

    std::optional<real_solution_t> getRealSolutions() const {
        auto delta = equation.getDelta();
        if (delta >= 0) {
            double root = std::sqrt(delta);
            double solution_1 = (-equation.getB() + root) / 2 / equation.getA();
            double solution_2 = (-equation.getB() - root) / 2 / equation.getA();
            return std::make_pair(solution_1, solution_2);
        }
        return std::nullopt;
    }

    std::optional<complex_solution_t> getComplexSolutions() const {
        auto delta = equation.getDelta();
        if (delta < 0) {
            double root = sqrt(-delta);
            return std::make_pair(
                    std::complex(-equation.getB() / 2 / equation.getA(), -root / 2 / equation.getA()),
                    std::complex(-equation.getB() / 2 / equation.getA(), +root / 2 / equation.getA()));

        }
        return std::nullopt;
    }

    void printSolutions() const {
        if (auto realSolutions = getRealSolutions()) {
            std::cout << "Equation has real solutions:\n"
                      << "x1 = " << realSolutions->first << "\tx2 =  " << realSolutions->second;
        } else if (auto complexSolutions = getComplexSolutions()) {
            std::cout << "Equation has complex solutions:\n"
                      << "z1 = " << complexSolutions->first << "\tz2 = " << complexSolutions->second;
        }
    }

};



int main() {
    Equation equation(21, 2, 1);
    std::cout << equation;
    EquationSolver equationSolver(equation);
    equationSolver.printSolutions();

    equation.setCoefficients(1,-2,1);
    std::cout << '\n' << equation;
    equationSolver.printSolutions();

    return 0;
}
