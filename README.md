# iSEL-Open-Source: Alkaline Water Electrolysis System Simulation Python Module

신재생에너지 연계 P2G 공정의 디지털 트윈: 수소생산을 위한 알칼라인 수전해 시스템

저자: 신용범, 오종연, 서혜원, 신동일

Abstract
신재생에너지는 변동성이 크고 간헐적이므로 전력망의 안정성 유지를 위한 장기적 대용량 전기 저장이 필요하며, 이를 위해 Power-to-Gas (P2G) 기술이 활용된다. P2G 기술을 통한 수소 생산 공정 및 저장장치의 위험성 관리 및 평가는 필수적이며, 신재생에너지 도입율이 증가하는 경우 더 많은 초과 에너지를 효율적으로 저장, 공급하기 위한 안정적이고 안전한 시스템 운영이 필요하다. 본 연구는 새만금 재생에너지 국가종합실증연구단지에 구축될 500 kW급 알칼라인 수전해시스템의 디지털 트윈 구축을 통해 고도화된 시스템 운전 및 관리 기술 개발을 목표로 하고 있으며, first-principles based 모델을 토대로 실제 운전 데이터와 결합하여 모델 parameter를 업데이트하는 디지털 트윈을 개발하고자 한다. 디지털 트윈을 구성하는 핵심 모델은 다음과 같다: 1) 온도와 압력에 따른 reversible voltage를 예측하는 열역학 모델, 2) 온도, 압력, 전류, 재료 등에 따른 셀 전압을 예측하는 전기화학 모델, 3)수전해시스템의 온도 변화 예측 모델, 4) 수소 생산량 및 산소내 수소 농도를 예측하는 성능 예측 모델. 디지털 트윈 구축을 위해 열역학 모델 4개와 전기화학 모델 3개를 기반으로 Python 모듈로 구축하였으며, 추후 온도 예측 모델과 성능 예측 모델을 Python 기반 모듈로 추가 구축하여 GitHub에 오픈소스로 공개하고, 이를 실제 운전 데이터와 결합하여 500 kW급 알칼라인 수전해시스템의 디지털 트윈을 구축할 예정이다.

- Reversible Voltagge
![image](https://user-images.githubusercontent.com/82799346/115238714-f71d0400-a158-11eb-8557-8b66dbd09492.png)

- Cell Voltage
![image](https://user-images.githubusercontent.com/82799346/115238748-fdab7b80-a158-11eb-863d-8166717d896a.png)

![image](https://user-images.githubusercontent.com/82799346/115239545-d4d7b600-a159-11eb-9523-1ba296898bc0.png)

- Thermal Model - 구축 예정
- Hydrogen Production Model - 구축 예정

알칼라인 수전해시스템의 DT 구축 위한 모델 확보: 
열역학 모델 4개와 전기화학 모델 3개를 수집하여 그림 6의 Python 모듈로 구축
대상 수전해시스템에 적합한 모델 선택 통한 customization 제공.
추후 온도 예측 모델 및 성능 예측 모델을 Python 모듈로 추가하여 전체 시스템에 대한 동적 시뮬레이션 모듈을 구축하고, 실제 데이터와 융합된 디지털 트윈 구축 예정.
