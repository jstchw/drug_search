import axios, {AxiosResponse} from 'axios'
import React from "react"
import { Nav } from "react-bootstrap"
import { animated, useTransition } from "@react-spring/web"


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

    const transitions = useTransition(selectedFeature, {
        from: { opacity: 0, transform: 'translate3d(100%,0,0)' },
        enter: { opacity: 1, transform: 'translate3d(0%,0,0)' },
        leave: { opacity: 0, transform: 'translate3d(-50%,0,0)' },
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

    React.useEffect(() => {
        console.log(selectedFeature)
    }, [selectedFeature])

    return (
        <React.Fragment>
            <Nav variant={'underline'} className={'d-flex justify-content-center mb-4'}>
                {uniqueFeatures && uniqueFeatures.map((feature: string, index: number) =>
                    <Nav.Link key={index} onClick={() => setSelectedFeature(feature)}
                              className={selectedFeature === feature ? 'active' : ''}>{feature}</Nav.Link>
                )}
            </Nav>

            {transitions((style, item) => item && (
                <animated.div style={style}>
                    {causalInfo && causalInfo.map((item: CausalInfoType, index: number) => {
                        if (item.feature.toUpperCase() === selectedFeature) {
                            return (
                                <div key={index}>
                                    <h1>Value: {item.value}</h1>
                                    <h1>Z-score: {item.z_score}</h1>
                                </div>
                            )
                        }
                        return null
                    })}
                </animated.div>
            ))}
        </React.Fragment>
    )
}

export default CausalInterface