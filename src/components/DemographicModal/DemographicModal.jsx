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

    for (let key in data) {
        // Remove unnecessary properties
        delete data[key].def;
        delete data[key].age;

        let termCountDict = data[key].termCountDict;
        // Convert the object to an array of entries and sort it
        let sortedDict = Object.entries(termCountDict).sort((a, b) => b[1] - a[1]);
        // Keep only the top 10 items and convert them back into an object
        let popularTerms = {};
        for (let [term, count] of sortedDict.slice(0, maxTerms)) {
            popularTerms[term] = count;
        }

        // Replace the original termCountDict with the sorted one
        data[key].termCountDict = popularTerms;
    }
    return data
}

const DemographicModalInfo = ({word}) => {
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
    }, [word]);

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
                    <div className={'fs-6 text-center'}><strong>{count.toLocaleString()}</strong> cases recorded</div>
                </Popover.Body>
            </Popover>
        )

        return (
            <OverlayTrigger trigger={['hover', 'focus']} placement={'top'} overlay={popover}>
                <Badge className={'mx-1 no-pointer-badge'}>{name}</Badge>
            </OverlayTrigger>
        )
    }


    return (
        <React.Fragment>
            {drugInfo && drugInfo.length > 0 &&
                <Row>
                    <Col xs={4}>
                        <DrugDescription
                            drugInfo={drugInfo[0]}
                            style={{border: "none"}}
                        />
                    </Col>
                    <Col xs={4}>
                        <DrugAccordion drugInfo={drugInfo[0]}/>
                    </Col>
                    <Col xs={4}>
                        <h4><Badge>AI</Badge>&nbsp;Demographic Summary</h4>
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
                </Row>
            }
            <Col className={'text-center'}>
                {isLoading ?
                    <Spinner animation="border" role="status">
                        <span className="visually-hidden">Loading...</span>
                    </Spinner>
                    :
                    (demographicInfo && Object.keys(demographicInfo).length > 0) ?
                        chunkArray(Object.entries(demographicInfo), 2).map((chunk, index) => (
                            <Row key={index}>
                                {chunk.map(([key, value], subIndex) => (
                                    <Col className={'text-center'} key={subIndex}>
                                        <h4>{getDemographicDef(key, value.totalCount, value.def)}</h4>
                                        <ApexChart
                                            eventDict={value.termCountDict}
                                            totalCount={value.totalCount}
                                            type={'searched_group'} />
                                    </Col>
                                ))}
                            </Row>
                        ))
                        :
                        <p>No demographic info available</p>
                }
            </Col>
        </React.Fragment>
    )
}

const DemographicModal = ({ show, handleClose, selectedWord}) => {
    const [word, setWord] = useState(null)

    // Make the word not all caps and check if the word has changed
    useEffect(() => {
        if (selectedWord) {
            const newWord = selectedWord.charAt(0).toUpperCase() + selectedWord.slice(1).toLowerCase();
            setWord(prevWord => newWord !== prevWord ? newWord : prevWord);
        }
    }, [selectedWord]);

    const handleOnHide = () => {
        handleClose();
        setWord(null);
    }


    return (
    <Modal centered show={show} onHide={handleOnHide} size={"xl"}>
      <Modal.Header closeButton>
          {word && <Modal.Title>{word}</Modal.Title>}
      </Modal.Header>
      <Modal.Body>
            {word && <DemographicModalInfo word={word}/>}
      </Modal.Body>
    </Modal>
    )
}

export default DemographicModal