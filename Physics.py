
class Solver1:
    """Это класс отвечающий за расчет коэффициентов учета поперечного краевого эффекта.
    Пока в тестовом режиме он называется Boltoncalculate. Он расчитывает только коэффициент Болтона
    """


    def __init__(self, data):
        self.a = data # Массив с описанием данных
        self.L = 1

    def calculate_resistance(self):
        """Рассчитываем взамное магнитное сопротивление вевей
        """
        size_j = len(self.a[0])
        size_i = len(self.a)
        size = (size_i-1)*(size_j-1)
        rmn_left = [0 for j in range(size)]
        rmn_right = [0 for j in range(size)]
        rmt_up = [0 for j in range(size)]
        rmt_down = [0 for j in range(size)]
        rmn = [0 for j in range(size)]
        rmt = [0 for j in range(size)]
        counter = 0
        for i in range(size_i-1):
            for j in range(size_j-1):
                rmn_left[counter] = (self.a[i][j+1].center.y - self.a[i][j].center.y) / \
                                    (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu() + 1 / self.a[i][j+1].mu())

                rmn_right[counter] = (self.a[i+1][j+1].center.y - self.a[i+1][j].center.y) / \
                                     (4 * self.a[i+1][j].width() * self.L) * (1 / self.a[i+1][j].mu() + 1
                                                                              / self.a[i+1][j+1].mu())

                rmt_down[counter] = (self.a[i+1][j].center.x - self.a[i][j].center.x) / \
                                    (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu() + 1
                                                                            / self.a[i][j+1].mu())

                rmt_up[counter] = (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) / \
                                  (4 * self.a[i][j + 1].height() * self.L) * \
                                  (1 / self.a[i][j + 1].mu() + 1 / self.a[i][j].mu())

                rmn[counter] = ((self.a[i][j + 1].center.y - self.a[i][j].center.y) /
                               (4 * self.a[i][j].width() * self.L) * (1 / self.a[i][j].mu() + 1 / self.a[i][j+1].mu())
                               + (self.a[i + 1][j + 1].center.y - self.a[i + 1][j].center.y) /
                               (4 * self.a[i + 1][j].width() * self.L) * (1 / self.a[i + 1][j].mu() + 1 /
                                                                          self.a[i+1][j+ 1].mu()))

                rmt[counter] = ((self.a[i + 1][j].center.x - self.a[i][j].center.x) /
                               (4 * self.a[i][j].height() * self.L) * (1 / self.a[i][j].mu()
                               + 1 / self.a[i][j + 1].mu())
                               + (self.a[i + 1][j + 1].center.x - self.a[i][j + 1].center.x) /
                               (4 * self.a[i][j + 1].height() * self.L) * (1 / self.a[i][j + 1].mu()
                                                                           + 1 / self.a[i][j].mu()))
                counter += 1
        return rmn_left, rmn_right, rmt_up, rmt_down, rmn, rmt





