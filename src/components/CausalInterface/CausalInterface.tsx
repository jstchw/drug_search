import axios, {AxiosResponse} from 'axios'
import React from "react"
import { Nav } from "react-bootstrap"
import { animated, useTransition } from "@react-spring/web"
import { Table } from "react-bootstrap"


interface CausalInfoType {
    feature: string,
    value: string,
    z_score: number,
}

const getCausalInfo = async (): Promise<CausalInfoType[] | null> => {
    return axios.get<CausalInfoType[]>('http://localhost:16000/api/get_causes', {responseType: 'json'})
        .then((response: AxiosResponse<CausalInfoType[]>) => {
            return response.data;
        })
        .catch((error) => {
            throw error;
        });
}

const CausalInterface = () => {
    const [causalInfo, setCausalInfo] = React.useState<CausalInfoType[] | null>(null)
    const uniqueFeatures = causalInfo && [...new Set(causalInfo.map((item: CausalInfoType) => item.feature.toUpperCase()))]
    const [selectedFeature, setSelectedFeature] = React.useState<string>('PSD')
    const [currentFeatureIndex, setCurrentFeatureIndex] = React.useState<number>(0)
    const prevFeatureIndex = React.useRef<number>(0)
    const direction = currentFeatureIndex > prevFeatureIndex.current ? 'right' : 'left'

    const transitions = useTransition(selectedFeature, {
        from: { opacity: 0, transform: direction === 'right' ? 'translate3d(100%,0,0)' : 'translate3d(-100%,0,0)' },
        enter: { opacity: 1, transform: 'translate3d(0%,0,0)' },
        exit: { opacity: 0, transform: direction === 'right' ? 'translate3d(-100%,0,0)' : 'translate3d(100%,0,0)' },
    })

    React.useEffect(() => {
        const fetchData = (): void  => {
            getCausalInfo()
                .then((data: CausalInfoType[] | null) => {
                    setCausalInfo(data)
                })
                .catch((error) => {
                    console.log(error)
                })
        }
        fetchData()
    }, [])

    const handleFeatureChange = (index: number, feature: string): void => {
        prevFeatureIndex.current = currentFeatureIndex
        setCurrentFeatureIndex(index)
        setSelectedFeature(feature)
    }

    return (
        <React.Fragment>
            <Nav variant={'underline'} className={'d-flex justify-content-center mb-4'}>
                {uniqueFeatures && uniqueFeatures.map((feature: string, index: number) =>
                    <Nav.Link key={index} onClick={() => handleFeatureChange(index, feature)}
                              className={selectedFeature === feature ? 'active' : ''}>{feature}</Nav.Link>
                )}
            </Nav>

            {transitions((style, item) => item && (
                <animated.div style={style}>
                    <Table striped bordered hover>
                        <thead>
                            <tr>
                                <th>Value</th>
                                <th>Z-score</th>
                            </tr>
                        </thead>
                        {causalInfo && causalInfo.map((item: CausalInfoType, index: number) => {
                            if (item.feature.toUpperCase() === selectedFeature && item.z_score >= 7) {
                                return (
                                    <tbody key={index}>
                                        <tr>
                                            <td>{item.value}</td>
                                            <td>{Math.round(item.z_score * 100) / 100}</td>
                                        </tr>
                                    </tbody>
                                )
                            }
                            return null
                        })}
                    </Table>
                </animated.div>
            ))}
        </React.Fragment>
    )
}

export default CausalInterface