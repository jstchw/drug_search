import React from 'react'
import {Badge, Col, Modal, OverlayTrigger, Placeholder, Popover, Row, Spinner} from 'react-bootstrap'
import {useEffect, useState} from "react";
import DrugDescription from "../DrugDescription/DrugDescription";
import DrugAccordion from "../DrugAccordion/DrugAccordion";
import useDrugInfo from "../../hooks/useDrugInfo";
import ApexChart from "../ApexChart/ApexChart";
import { getSideEffectsForDemographics } from "../../services/FDA_Request";
import { searchTypes } from "../OptionModal/OptionModal";

// Function to split the array in chunks,
const chunkArray = (array, size) => {
    const chunked_arr = [];
    let index = 0;
    while (index < array.length) {
        chunked_arr.push(array.slice(index, size + index));
        index += size;
    }
    return chunked_arr;
};

const prepareDataForGPT = (data) => {
    let maxTerms = 10;
    let gptData = JSON.parse(JSON.stringify(data))

    for (let key in gptData) {
        if(gptData.hasOwnProperty(key)) {
            delete gptData[key].def
            let sortedTerms = Object.entries(gptData[key].termCountDict).sort((a, b) => b[1] - a[1])
            gptData[key].termCountDict = sortedTerms.slice(0, maxTerms).reduce((obj, [k, v]) => ({...obj, [k]: v}), {})
        }
    }
    return gptData
}

const DemographicModalInfo = ({word, isMobile}) => {
    const wordRef = React.useRef(null)
    const drugInfo = useDrugInfo([word], searchTypes[0].value)
    const [demographicInfo, setDemographicInfo] = useState([])
    const [demographicSummary, setDemographicSummary] = useState(null)
    const [isLoading, setIsLoading] = useState(false)
    let placeholderSizes = [6, 4, 4, 5, 3, 4, 4, 5, 3, 2, 4, 6]
    placeholderSizes = placeholderSizes.concat(placeholderSizes)

    // Get demographic info only if the word has changed
    useEffect(() => {
        if (word !== wordRef.current) {
            wordRef.current = word;
            setIsLoading(true);
            getSideEffectsForDemographics(word)
                .then((res) => {
                    setDemographicInfo(res);
                })
                .catch((err) => {
                    console.log(err);
                })
                .finally(() => {
                    setIsLoading(false);
                });
        }
    }, [word])

    useEffect(() => {
        console.log(demographicInfo)
    })

    // Get response from the GPT model if the demographic info has changed
    useEffect(() => {
        if(Object.keys(demographicInfo).length > 0) {
            const cachedResponse = localStorage.getItem(JSON.stringify(demographicInfo));
            if (cachedResponse) {
                setDemographicSummary(JSON.parse(cachedResponse));
            } else {
                fetch('api/get_summary', {
                    method: 'POST',
                    headers: {
                        'Content-Type': 'application/json'
                    },
                    body: JSON.stringify(prepareDataForGPT(demographicInfo))
                })
                    .then(res => res.json())
                    .then(data => {
                        setDemographicSummary(data)
                        localStorage.setItem(JSON.stringify(demographicInfo), JSON.stringify(data))
                    })
                    .catch(err => console.log(err))
            }
        }
    }, [demographicInfo])


    const getDemographicDef = (name, count, def) => {
        const popover = (
            <Popover>
                <Popover.Body>
                    <div className={'fs-6 text-center'}>{def}</div>
                    <div className={'fs-6 text-center'}><strong>{count.toLocaleString()}</strong> cases recorded
                    </div>
                </Popover.Body>
            </Popover>
        )

        return (
            <OverlayTrigger trigger={['hover', 'focus']} placement={'top'} overlay={popover}>
                <Badge className={'mx-1 no-pointer-badge'}>{name}</Badge>
            </OverlayTrigger>
        )
    }

    const showAIWarning = () => {
        const popover = (
            <Popover>
                <Popover.Body>
                    <div className={'fs-6 text-center'}>
                        This content is generated using GPT 3.5 model and may contain false information.
                        Always cross-verify with reliable sources and consult professionals for important decisions.
                    </div>
                </Popover.Body>
            </Popover>
        )

        return (
            <OverlayTrigger trigger={['hover', 'focus']} placement={isMobile ? 'right' : 'left'} overlay={popover}>
                <Badge className={'mx-1 no-pointer-badge'}>AI</Badge>
            </OverlayTrigger>
        )
    }

    const DrugSummaryObject = ({drugInfo, demographicSummary, placeholderSizes, width}) => {
        return (
            <React.Fragment>
                {drugInfo && drugInfo.length > 0 && (
                    <React.Fragment>
                        <Col xs={width}>
                            <DrugDescription
                                drugInfo={drugInfo[0]}
                                style={{border: "none"}}
                            />
                        </Col>
                        <Col xs={width} className={isMobile ? 'mb-4' : ''}>
                            <DrugAccordion drugInfo={drugInfo[0]}/>
                        </Col>
                    </React.Fragment>
                )}
                <Col xs={width}>
                    <h4>{showAIWarning()}Demographic Summary</h4>
                    {demographicSummary ? (
                        <p>{demographicSummary.choices[0].message.content}</p>
                    ) : (
                        <Placeholder as="p" animation={'glow'}>
                            {placeholderSizes.map((size, index) => (
                                <React.Fragment key={index}>
                                    <Placeholder xs={size} />{' '}
                                </React.Fragment>
                            ))}
                        </Placeholder>
                    )}
                </Col>
            </React.Fragment>
        )
    }

    return (
        <React.Fragment>
            {!isMobile ?
                <Row className={!drugInfo.length > 0 ? 'd-flex justify-content-center text-center' : ''}>
                    <DrugSummaryObject drugInfo={drugInfo} demographicSummary={demographicSummary}
                                       placeholderSizes={placeholderSizes} width={4}/>
                </Row>
                :
                <DrugSummaryObject drugInfo={drugInfo} demographicSummary={demographicSummary}
                                   placeholderSizes={placeholderSizes} width={12}/>
            }

            <Col className={'text-center'}>
                {isLoading ?
                    <Spinner animation="border" role="status">
                        <span className="visually-hidden">Loading...</span>
                    </Spinner>
                    :
                    (demographicInfo && Object.keys(demographicInfo).length > 0) ?
                        chunkArray(Object.entries(demographicInfo), 2).map((chunk, index) => {
                            const Wrapper = isMobile ? React.Fragment : Row
                            return (
                                <Wrapper key={index}>
                                    {chunk.map(([key, value], subIndex) => (
                                        <Col className={'text-center'} key={subIndex}>
                                            <h4>{getDemographicDef(key, value.totalCount, value.def)}</h4>
                                            <ApexChart
                                                eventDict={value.termCountDict}
                                                totalCount={value.totalCount}
                                                type={'searched_group'}/>
                                        </Col>
                                    ))}
                                </Wrapper>
                            )
                        })
                        :
                        <p>No demographic info available</p>
                }
            </Col>
        </React.Fragment>
    )
}

const DemographicModal = ({ show, handleClose, selectedWord, isMobile}) => {
    const [word, setWord] = useState(null)
    const [shouldRenderContent, setShouldRenderContent] = useState(false)

    // Make the word not all caps and check if the word has changed
    useEffect(() => {
        if (selectedWord) {
            const newWord = selectedWord.charAt(0).toUpperCase() + selectedWord.slice(1).toLowerCase()
            setWord(newWord)
        }
    }, [selectedWord])

    useEffect(() => {
        if (show) {
            setShouldRenderContent(true);
        } else {
            const timer = setTimeout(() => {
                setShouldRenderContent(false);
            }, 300); // Delay for closing animation
            return () => clearTimeout(timer);
        }
    }, [show]);


    const handleOnHide = () => {
        handleClose()
    }


    return (
    <Modal centered show={show} onHide={handleOnHide} size={"xl"}>
      <Modal.Header closeButton>
          {word && <Modal.Title>{word}</Modal.Title>}
      </Modal.Header>
      <Modal.Body>
            {shouldRenderContent && word && <DemographicModalInfo word={word} isMobile={isMobile}/>}
      </Modal.Body>
    </Modal>
    )
}

export default DemographicModal