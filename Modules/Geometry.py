from Modules.DiscretizationNew import Coord

class Geometry:
    """Это класс предназначен для описания основных улов геометрии"""


    gNodes = {} # инициализируем словарь для узлов геометрии


    def __init__(self, label=''):
        self.label = label  # определяем имя геометрии


    def rect(self, height, width, intial=('intialX', 'intialY'), label='rect'):
        """Создаем прямоугольник с высотой height шириной width
        левый нижний угол intial = (x,y)
        label - имя узла
        """
        start = Coord(intial[0], intial[1])
        distr = [[start, Coord(start.x+width, start.y+height)]]
        self.gNodes[label] = distr

        return self


    def coils(self, height, width, intial=('intialX', 'intialY'), step=('stepX', 'stepY'), amount=1, label='coil'):
        """Функция создает прямоугольную катушку
        Аргументы:
        height - высота катушки по y
        width - ширина катушки по x
        intialX, initialY - начальная координата для прямоугольника (нижний левый угол)
        """
        distr = []
        start = Coord(intial[0], intial[1])
        step = Coord(step[0], step[1])
        for i in range(amount):
            x0 = start.x + i * step.x
            x1 = start.x + i * step.x + width
            y0 = start.y + i * step.y
            y1 = start.y + i * step.y + height
            distr.append([Coord(x0, y0), Coord(x1, y1)])
        self.gNodes[label] = distr
        return self

    def copy(self, object, gNodesLabel='test', step=('xStep', 'yStep'), number=1, label='copy'):
        """Копирует указаный узел геометрии с помощб лейбла (gNodesLabel)
        Копирует с шагом step =(x-step, y-step)
        number - сколько раз копировать
        label - новый лейбл скопированного узла
        """
        distr = []
        step = Coord(step[0], step[1])
        for i in range(number):
            for j in object.gNodes[gNodesLabel]:
                xNew1 = j[0].x + i * step.x
                yNew1 = j[0].y + i * step.y
                xNew2 = j[1].x + i * step.x
                yNew2 = j[1].y + i * step.y

                distr.append([Coord(xNew1, yNew1),Coord(xNew2, yNew2)])
        self.gNodes[label] = distr
        return self

    def delete(self, label=''):
        """Удаляет узел геометрии по указаному лейблу"""
        self.gNodes.pop(label)

    def setParamInductor(self, yokeHeight=0, slotHeight=0, toothWidth=0, slotWidth=0, slotNumber=0, label=''):
        setParamLocal = {}
        setParamLocal['yokeHeight'] = yokeHeight
        setParamLocal['slotHeight'] = slotHeight
        setParamLocal['toothWidth'] = toothWidth
        setParamLocal['slotWidth'] = slotWidth
        setParamLocal['slotNumber'] = slotNumber
        self.setParamNode= {label : setParamLocal}

    def inductor(self, *param, intialCoord=('X', 'Y'), labelSetings=None, label='inductor'):
        """Функция создает геометрию индуктора, ярмо и зубцы
        intialCoord - координаты нижнего левого угла ярма
        label - название элемента
        """
        if labelSetings == None:
            yokeHeight = param[0]
            slotHeight = param[1]
            slotWidth = param[2]
            toothWidth = param[3]
            slotNumber = param[4]
        else:
            yokeHeight = self.setParamNode[labelSetings]['yokeHeight']
            slotHeight = self.setParamNode[labelSetings]['slotHeight']
            slotWidth = self.setParamNode[labelSetings]['toothWidth']
            toothWidth = self.setParamNode[labelSetings]['slotWidth']
            slotNumber = self.setParamNode[labelSetings]['slotNumber']

        distr = []
        start = Coord(intialCoord[0], intialCoord[1])
        step = toothWidth + slotWidth
        #print(step)
        yokeWidth = step * slotNumber + toothWidth
        yoke = [start, Coord(start.x + yokeWidth, start.y + yokeHeight)]
        distr.append(yoke)
        for i in range(slotNumber+1):
            x0 = start.x + i * step
            x1 = start.x + slotHeight + i * step
            y0 = start.y + yokeHeight
            y1 = start.y + yokeHeight + slotHeight
            distr.append([Coord(x0, y0), Coord(x1, y1)])

        self.gNodes[label] = distr
        return self

