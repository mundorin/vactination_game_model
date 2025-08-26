import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve
from models import SIR_V

class SIR_V_CHOICE:
    def __init__(self, s, i, r, r_coef, c_i, c_v, treshhold=0.0001, smoothness=0.001, max_iter = 100000):
        if s + i + r != 1:
            print("error s + i + r != 1")
        self.s = s
        self.i = i
        self.r = r
        self.b = r_coef * smoothness
        self.y = smoothness
        self.c_i = c_i
        self.c_v = c_v

        self.treshhold = treshhold
        self.r_coef = r_coef
        self.smoothness = smoothness
        self.max_iter = max_iter

        self.time_step = 1
        self.timestamp = 0
        if self.i == 0:
            self.end = True
        else:
            self.end = False

        self.v, self.p_infect, self.score, _, _ = self.get_nash_equilibrium()
        self.v_optimal, self.p_infect_optimal, self.score = self.get_pareto_equilibrium()
        self.s -= self.v

        self.model_info = self.make_simulation()

    def make_simulation(self):
        max_iter = self.max_iter
        param = {
            's': self.s,
            'i': self.i,
            'r': self.r,
            'v': self.v,
            'b': self.b,
            'y': self.y,
            'n': 0,
            'treshhold': self.treshhold,
        }
        model = SIR_V(*param.values())

        row_info = model.get_info()
        model_info = pd.DataFrame({
            "timestamp": [row_info['timestamp']],
            "s": [row_info["s"]],
            "i": [row_info["i"]],
            "r": [row_info["r"]],
            "v": [row_info["v"]], })

        while not row_info["end"]:
            if row_info["timestamp"] > max_iter:
                print("stopped due to exceeding the maximum number of iterations")
                break
            model.iteretion()
            row_info = model.get_info()
            model_info.loc[len(model_info)] = {
                "timestamp": row_info['timestamp'],
                "s": row_info["s"],
                "i": row_info["i"],
                "r": row_info["r"],
                "v": row_info["v"], }

        return model_info

    def diplay_plots(self):
        plt.plot(self.model_info.s)
        plt.plot(self.model_info.i)
        plt.plot(self.model_info.r)
        plt.plot(self.model_info.v, color = "blue")
        plt.show()

    def iteretion(self):
        self.timestamp += self.time_step

        recovereds = self.y * self.i

        self.i -= recovereds
        self.r += recovereds

        sick = self.b * self.s * self.i

        self.s -= sick
        self.i += sick
        if not self.end and self.i <= self.treshhold:
            self.end = True

    def get_info(self):
        info = {
            "s": self.s,
            "i": self.i,
            "r": self.r,
            "r_coef": self.r_coef,
            "v": self.v,
            "n": self.n,
            "timestamp": self.timestamp,
            "end": self.end,
            "smoothness": self.smoothness,
            "p_infect": self.p_infect,
        }
        return info

    def get_social_score(self, v, r):
        res = v * self.c_v + r * self.c_i
        return res

    def get_infection_prob(self, v):
        """
        Вероятность заражения как функция от уровня вакцинации
        Решает уравнение: (1 - r_inf - v) = (1 - v) * exp(-R0 * r_inf)
        """
        R0 = self.r_coef

        # Порог коллективного иммунитета
        herd_threshold = 1 - 1 / R0

        # Если достигнут коллективный иммунитет
        if v >= herd_threshold:
            return 0.0

        def equation(r_inf):
            return (1 - r_inf - v) - (1 - v) * np.exp(-R0 * r_inf)

        # Начальное приближение - улучшенная версия
        initial_guess = min(0.8 * (1 - v), 0.5)  # Более устойчивое начальное приближение

        try:
            # Используем более надежный метод
            r_inf_solution = fsolve(equation, initial_guess, full_output=True, xtol=1e-6)

            if r_inf_solution[2] == 1:  # successful solution
                result = r_inf_solution[0][0]
            else:
                # Если fsolve не сходится, используем бинарный поиск
                result = self._solve_infection_prob_binary(v, R0)
        except:
            result = self._solve_infection_prob_binary(v, R0)

        # Обеспечиваем, что решение в пределах [0, 1-v]
        return max(0, min(result, 1 - v))

    def _solve_infection_prob_binary(self, v, R0, tolerance=1e-6):
        """
        Решение уравнения бинарным поиском (более надежно)
        """

        def equation(r_inf):
            return (1 - r_inf - v) - (1 - v) * np.exp(-R0 * r_inf)

        left, right = 0.0, 1 - v

        # Проверяем границы
        if equation(left) * equation(right) > 0:
            # Если на границах одинаковый знак, выбираем приближение
            return 0.5 * (1 - v) if abs(equation(left)) < abs(equation(right)) else 0.8 * (1 - v)

        # Бинарный поиск
        for _ in range(100):
            mid = (left + right) / 2
            f_mid = equation(mid)

            if abs(f_mid) < tolerance:
                return mid

            if f_mid > 0:
                left = mid
            else:
                right = mid

        return (left + right) / 2

    def get_pareto_equilibrium(self):
        """
        Находит Парето-оптимальное равновесие с улучшенной стабильностью
        """

        def objective(v):
            r = self.get_infection_prob(v)
            return self.get_social_score(v, r)

        # Сеточный поиск с меньшим количеством точек для скорости
        v_values = np.linspace(0, 1, 50)
        best_score = float('inf')
        best_v = 0.0

        for v in v_values:
            try:
                score = objective(v)
                if score < best_score:
                    best_score = score
                    best_v = v
            except:
                continue

        # Уточняем вокруг найденного минимума
        try:
            refine_range = 0.2
            v_min = max(0, best_v - refine_range)
            v_max = min(1, best_v + refine_range)

            result = minimize_scalar(
                objective,
                bounds=(v_min, v_max),
                method='bounded',
                options={'xatol': 1e-4, 'maxiter': 50}
            )
            if result.fun < best_score:
                best_v, best_score = result.x, result.fun
        except:
            pass

        optimal_r = self.get_infection_prob(best_v)
        return best_v, optimal_r, best_score

    def get_nash_equilibrium(self):
        """
        Находит равновесие Нэша (где индивидуально оптимально вакцинироваться)
        Условие: r * c_i = c_v
        Возвращает: (nash_v, nash_r, nash_condition)
        """
        # Решаем уравнение: get_infection_prob(v) * c_i = c_v
        # => base_infection_rate * (1 - v) * c_i = c_v
        # => 1 - v = c_v / (base_infection_rate * c_i)
        # => v = 1 - c_v / (base_infection_rate * c_i)

        # Проверяем существование решения
        critical_value = self.c_v / (self.c_i)

        if critical_value < 0:
            # c_v < 0 (невозможно) или c_i < 0 (невозможно)
            nash_v = 1.0  # все вакцинируются
            nash_condition = "invalid costs"

        elif critical_value <= 1:
            # Внутреннее равновесие Нэша
            nash_v = max(0, min(1, 1 - critical_value))
            nash_condition = "mixed equilibrium"

        else:  # critical_value > 1
            # c_v слишком высока или c_i слишком низка
            nash_v = 0.0  # никто не вакцинируется
            nash_condition = "no vaccination equilibrium"

        nash_r = self.get_infection_prob(nash_v)
        actual_condition = nash_r * self.c_i
        score = self.get_social_score(nash_v, nash_r)

        return nash_v, nash_r, score, nash_condition, actual_condition

    def analyze_equilibria(self):
        """
        Анализ равновесий Нэша и Парето с строгим выводом
        """
        # Вычисление равновесий
        pareto_v, pareto_r, pareto_score = self.get_pareto_equilibrium()
        nash_v, nash_r, nash_score, nash_condition, actual_condition = self.get_nash_equilibrium()

        # Анализ эффективности
        efficiency_gap = nash_score - pareto_score
        is_efficient = abs(efficiency_gap) < 1e-6

        # Преобразование numpy типов в float и округление
        def convert_and_round(value):
            if hasattr(value, 'item'):  # для numpy типов
                return round(float(value.item()), 4)
            return round(float(value), 4)

        nash_v = convert_and_round(nash_v)
        nash_r = convert_and_round(nash_r)
        nash_score = convert_and_round(nash_score)
        pareto_v = convert_and_round(pareto_v)
        pareto_r = convert_and_round(pareto_r)
        pareto_score = convert_and_round(pareto_score)
        efficiency_gap = convert_and_round(efficiency_gap)
        actual_condition_float = convert_and_round(actual_condition) if hasattr(actual_condition, 'item') else round(
            float(actual_condition), 4)

        # Строгий вывод информации
        print("=" * 80)
        print("АНАЛИЗ РАВНОВЕСИЙ МОДЕЛИ ВАКЦИНАЦИИ")
        print("=" * 80)

        print(f"\nПАРАМЕТРЫ МОДЕЛИ:")
        print(f"R₀: {self.r_coef:.2f}")
        print(f"c_v: {self.c_v:.1f}")
        print(f"c_i: {self.c_i:.1f}")
        print(f"Порог коллективного иммунитета: {1 - 1 / self.r_coef:.3f}")

        print(f"\nРАВНОВЕСИЕ НЭША:")
        print(f"v = {nash_v:.4f}")
        print(f"r = {nash_r:.4f}")
        print(f"Издержки = {nash_score:.4f}")
        print(f"Тип: {nash_condition}")
        print(f"Проверка условия: {actual_condition_float:.3f} vs {self.c_v:.1f}")

        print(f"\nПАРЕТО-ОПТИМУМ:")
        print(f"v = {pareto_v:.4f}")
        print(f"r = {pareto_r:.4f}")
        print(f"Издержки = {pareto_score:.4f}")

        print(f"\nСРАВНИТЕЛЬНЫЙ АНАЛИЗ:")
        print(f"Разница издержек: {efficiency_gap:.4f}")
        efficiency_status = "Эффективно" if is_efficient else "Неэффективно"
        print(f"Эффективность: {efficiency_status}")

        if not is_efficient:
            vaccination_gap = convert_and_round(abs(nash_v - pareto_v))
            if nash_v < pareto_v:
                print(f"Причина неэффективности: недостаточная вакцинация")
                print(f"Дефицит вакцинации: {vaccination_gap:.4f}")
            else:
                print(f"Причина неэффективности: избыточная вакцинация")
                print(f"Избыток вакцинации: {vaccination_gap:.4f}")

        print("=" * 80)

        # Возврат структурированных данных с преобразованными типами
        return {
            'parameters': {
                'R0': float(self.r_coef),
                'c_v': float(self.c_v),
                'c_i': float(self.c_i),
                'herd_threshold': round(1 - 1 / float(self.r_coef), 4)
            },
            'nash_equilibrium': {
                'vaccination_rate': nash_v,
                'infection_prob': nash_r,
                'social_cost': nash_score,
                'equilibrium_type': nash_condition,
                'condition_check': f"{actual_condition_float:.3f} vs {float(self.c_v):.1f}"
            },
            'pareto_optimum': {
                'vaccination_rate': pareto_v,
                'infection_prob': pareto_r,
                'social_cost': pareto_score
            },
            'efficiency_analysis': {
                'cost_difference': efficiency_gap,
                'is_efficient': bool(is_efficient),
                'efficiency_status': efficiency_status,
                'vaccination_gap': vaccination_gap if not is_efficient else 0.0
            }
        }


class SIRV_CHOICE:
    def __init__(self, s, i, r, v, b, y, fi, n=10000, treshhold=0.0001):
        if s + i + r + v != 1:
            print("error s + i + r + v != 1")
        self.s = s
        self.i = i
        self.r = r
        self.v = v
        self.b = b
        self.y = y
        self.fi = fi
        self.time_step = 1
        self.timestamp = 0
        self.n = n
        self.treshhold = treshhold
        if self.i == 0:
            self.end = True
        else:
            self.end = False

    def iteretion(self):
        self.timestamp += self.time_step

        vactinatons = self.s * self.fi

        self.s -= vactinatons
        self.v += vactinatons

        recovereds = self.y * self.i

        self.i -= recovereds
        self.r += recovereds

        sick = self.b * self.s * self.i

        self.s -= sick
        self.i += sick

        if not self.end and self.i <= self.treshhold:
            self.end = True

    def get_info(self):
        info = {
            "s": self.s,
            "i": self.i,
            "r": self.r,
            "b": self.b,
            "y": self.y,
            "v": self.v,
            "fi": self.fi,
            "n": self.n,
            "timestamp": self.timestamp,
            "end": self.end,
        }
        return info
